package water.rapids;


import sun.misc.Unsafe;
import water.*;
import water.fvec.Chunk;
import water.fvec.Frame;
import water.fvec.NewChunk;
import water.fvec.Vec;
import water.nbhm.NonBlockingHashSet;
import water.nbhm.UtilUnsafe;
import water.util.Log;

import java.util.Arrays;
import java.util.HashMap;
import java.util.concurrent.atomic.AtomicInteger;


/**
 * GROUPBY: Single pass aggregation by columns.
 *
 * NA handling:
 *
 *  AGG.T_IG: case 0
 *    Count NA rows, but discard values in sums, mins, maxs
 *      FIRST/LAST return the first nonNA first/last, or NA if all NA
 *
 *  AGG.T_RM: case 1
 *    Count NA rows separately, discard values in sums, mins, maxs and compute aggregates less NA row counts
 *      FIRST/LAST treated as above
 *
 *  AGG.T_ALL: case 2
 *    Include NA in all aggregates -- any NA encountered forces aggregate to be NA.
 *      FIRST/LAST return first/last row regardless of NAs.
 *
 * Aggregates:
 *  MIN
 *  MAX
 *  MEAN
 *  COUNT
 *  SUM
 *  SD
 *  VAR
 *  COUNT_DISTINCT
 *  FIRST
 *  LAST
 *  Aggregations on time and numeric columns only.
 */
  public class ASTGroupBy extends ASTUniPrefixOp {
  // AST: (GB fr {cols} AGGS)
  //      (GB %k {#1;#3} (AGGS #2 "min" #4 "mean" #6))
  private long[] _gbCols; // group by columns
  private AGG[] _agg;
  ASTGroupBy() { super(null); }
  @Override protected String opStr() { return "GB"; }
  @Override protected ASTOp make() {return new ASTGroupBy();}
  protected ASTGroupBy parse_impl(Exec E) {
    AST ary = E.parse();
    if( ary instanceof ASTId ) ary = Env.staticLookup((ASTId)ary);

    // parse gby columns
    AST s=null;
    try {
      s=E.skipWS().parse();
      _gbCols=((ASTSeries)s).toArray();
      if(_gbCols.length > 1000 )
        throw new IllegalArgumentException("Too many columns selected. Please select < 1000 columns.");
    } catch (ClassCastException e) {
      assert s!=null;
      try {
        _gbCols = new long[]{(long)((ASTNum)s).dbl()};
      } catch (ClassCastException e2) {
        throw new IllegalArgumentException("Badly formed AST. Columns argument must be a ASTSeries or ASTNum");
      }
    }

    //parse AGGs
    _agg = ((AGG)E.parse())._aggs;

    ASTGroupBy res = (ASTGroupBy)clone();
    res._asts = new AST[]{ary};
    return res;
  }
  @Override protected void apply(Env e) {
    // only allow reductions on time and numeric columns
    Frame fr = e.popAry();

    GBTask p1 = new GBTask(_gbCols, _agg).doAll(fr);

    // build the output
    final int nGrps = p1._g.size();
    final int nCols = _gbCols.length+_agg.length;

    // dummy vec
    Vec v = Vec.makeZero(nGrps);

    // the names of columns
    String[] names = new String[nCols];
    String[][] domains = new String[nCols][];
    for( int i=0;i<_gbCols.length;++i) {
      names[i] = fr.name((int) _gbCols[i]);
      domains[i] = fr.domains()[(int)_gbCols[i]];
    }
    System.arraycopy(AGG.names(_agg),0,names, _gbCols.length,_agg.length);

    final G[] grps = p1._g.toArray(new G[nGrps]);
    final AGG[] agg=_agg;
    Frame f=new MRTask() {
      @Override public void map(Chunk[] c, NewChunk[] ncs) {
        int start=(int)c[0].start();
        for( int i=0;i<c[0]._len;++i) {
          G g = grps[i+start];
          int j=0;
          for(;j<g._ds.length;++j)
            ncs[j].addNum(g._ds[j]);

          for(int a=0; a<agg.length;++a) {
            byte type = agg[a]._type;
            switch( type ) {
              case AGG.T_N:  ncs[j++].addNum(g._N       );  break;
              case AGG.T_ND: ncs[j++].addNum(g._ND[a]   );  break;
              case AGG.T_F:  ncs[j++].addNum(g._f[a]    );  break;
              case AGG.T_L:  ncs[j++].addNum(g._l[a]    );  break;
              case AGG.T_MIN:ncs[j++].addNum(g._min[a]  );  break;
              case AGG.T_MAX:ncs[j++].addNum(g._max[a]  );  break;
              case AGG.T_AVG:ncs[j++].addNum(g._avs[a]  );  break;
              case AGG.T_VAR:ncs[j++].addNum(g._vars[a] );  break;
              case AGG.T_SD :ncs[j++].addNum(g._sdevs[a]);  break;
              case AGG.T_SUM:ncs[j++].addNum(g._sum[a]  );  break;
              case AGG.T_SS :ncs[j++].addNum(g._ss [a]  );  break;
              default:
                throw new IllegalArgumentException("Unsupported aggregation type: " + type);
            }
          }
        }
      }
    }.doAll(nCols,v).outputFrame(Key.make(),names,domains);
    p1._g=null;
    Keyed.remove(v._key);
    e.pushAry(f);
  }

  private static class GBTask extends MRTask<GBTask> {
    NonBlockingHashSet<G> _g;
    private long[] _gbCols;
    private AGG[] _agg;
    GBTask(long[] gbCols, AGG[] agg) { _gbCols=gbCols; _agg=agg; }
    @Override public void setupLocal() { _g = new NonBlockingHashSet<>(); }
    @Override public void map(Chunk[] c) {
      long start = c[0].start();
      byte[] naMethods = AGG.naMethods(_agg);
      for (int i=0;i<c[0]._len;++i) {
        G g = new G(i,c,_gbCols,_agg.length,naMethods);
        if( !_g.add(g) ) g=_g.get(g);
        // cas in COUNT
        long r=g._N;
        while(!G.CAS_N(g, r, r + 1))
          r=g._N;
        perRow(_agg,i,start,c,g);
      }
    }
    @Override public void reduce(GBTask t) {
      if( _g!=t._g ) {
        NonBlockingHashSet<G> l = _g;
        NonBlockingHashSet<G> r = t._g;
        if( l.size() < r.size() ) { l=r; r=_g; }  // larger on the left

        // loop over the smaller set of grps
        for( G rg:r ) {
          G lg = l.get(rg);
          if( !l.add(rg) ) {
            assert lg!=null;
            long R = lg._N;
            while (!G.CAS_N(lg, R, R + rg._N))
              R = lg._N;
            reduceGroup(_agg, lg, rg);
          }
        }
        _g=l;
        t._g=null;
      }
    }
    @Override public void closeLocal() {}
    @Override public void postGlobal() { H2O.submitTask(new ParallelPostGlobal(_g.toArray(new G[_g.size()]))).join(); }

    // task helper functions
    private static void perRow(AGG[] agg, int chkRow, long rowOffset, Chunk[] c, G g) { perRow(agg,chkRow,rowOffset,c,g,null); }
    private static void reduceGroup(AGG[] agg, G g, G that) { perRow(agg,-1,-1,null,g,that);}
    private static void perRow(AGG[] agg, int chkRow, long rowOffset, Chunk[] c, G g, G that) {
      byte type; int col;
      for (int i=0;i<agg.length;++i) {
        col = agg[i]._c;

        // update NA value for this (group, aggregate) pair:
        if( c!=null ) {
          if( c[col].isNA(chkRow) ) setNA(g,1L,i);
        } else {
          // reduce NA counts together...
          setNA(g,that._NA[i],i);
        }

        if( (type=agg[i]._type) == AGG.T_N ) continue; //immediate short circuit if COUNT

        // Do NA handling here for AGG.T_IG and AGG.T_RM
        if( c!= null )
          if( !agg[i].isAll() && c[col].isNA(chkRow) )
            continue;

        // build up a long[] of vals, to handle the case when c is and isn't null.
        // c is null in the reduce  of the MRTask
        long[] vals = new long[6]; // 6 cases in the switch, magic array for legibility in the switch
        long bits=-1;
        if( c!=null ) {
          if( c[col].isNA(chkRow) ) continue;
          bits = Double.doubleToRawLongBits(c[col].atd(chkRow));
        }
        vals[0] = c==null ? that._f[i] : chkRow+rowOffset;
        vals[1] = c==null ? that._l[i] : chkRow+rowOffset;
        vals[2] = c==null ? Double.doubleToRawLongBits(that._min[i]) : bits;
        vals[3] = c==null ? Double.doubleToRawLongBits(that._max[i]) : bits;
        vals[4] = c==null ? Double.doubleToRawLongBits(that._sum[i]) : bits;
        vals[5] = c==null ? Double.doubleToRawLongBits(that._ss[i] ) : bits;
        if( type == AGG.T_ND ) {
//          if( c==null ) g._nd._nd[i].addAll(that._nd._nd[i]);
//          else          g._nd._nd[i].add(c[col].atd(chkRow));
          continue;
        }

        switch( type ) {
          case AGG.T_F:   setFirst(g,vals[0],i);   break;
          case AGG.T_L:   setLast( g,vals[1],i);   break;
          case AGG.T_MIN: setMin(  g,vals[2],i);   break;
          case AGG.T_MAX: setMax(  g,vals[3],i);   break;
          case AGG.T_AVG: /* fall through */
          case AGG.T_SUM: setSum(  g,vals[4],i);   break;
          case AGG.T_VAR: /* fall through */
          case AGG.T_SD:
          case AGG.T_SS:  setSS(   g,vals[5],i);   break;
          default:
            throw new IllegalArgumentException("Unsupported aggregation type: " + type);
        }
      }
    }

    // all the CAS'ing helpers
    private static void setFirst(G g, long v, int c) {
      long o = g._f[c];
      while( v < o && !G.CAS_f(g,G.longRawIdx(c),o,v))
        o = g._f[c];
    }
    private static void setLast(G g, long v, int c) {
      long o = g._l[c];
      while( v > o && !G.CAS_l(g, G.longRawIdx(c), o, v))
        o = g._l[c];
    }
    private static void setMin(G g, long v, int c) {
      double o = g._min[c];
      double vv = Double.longBitsToDouble(v);
      while( vv < o && !G.CAS_min(g,G.doubleRawIdx(c),Double.doubleToRawLongBits(o),v))
        o = g._min[c];
    }
    private static void setMax(G g, long v, int c) {
      double o = g._max[c];
      double vv = Double.longBitsToDouble(v);
      while( vv > o && !G.CAS_max(g, G.doubleRawIdx(c), Double.doubleToRawLongBits(o), v))
        o = g._max[c];
    }
    private static void setSum(G g, long vv, int c) {
      double v = Double.longBitsToDouble(vv);
      double o = g._sum[c];
      while(!G.CAS_sum(g,G.doubleRawIdx(c),Double.doubleToRawLongBits(o),Double.doubleToRawLongBits(o+v)))
        o=g._sum[c];
    }
    private static void setSS(G g, long vv, int c) {
      double v = Double.longBitsToDouble(vv);
      double o = g._ss[c];
      while(!G.CAS_ss(g, G.doubleRawIdx(c), Double.doubleToRawLongBits(o), Double.doubleToRawLongBits(o + v * v)))
        o=g._ss[c];
    }
    private static void setNA(G g, long n, int c) {
      long o = g._NA[c];
      while(!G.CAS_NA(g,G.longRawIdx(c),o,o+n))
        o=g._NA[c];
    }

    // serialization
    @Override public AutoBuffer write_impl( AutoBuffer ab ) {
      ab.putA8(_gbCols);
      ab.put4(_agg.length);
      for (AGG a_agg : _agg) ab.put(a_agg);
      if( _g == null ) return ab.put4(0);
      ab.put4(_g.size());
      for( G g: _g) ab.put(g);
      return ab;
    }

    @Override public GBTask read_impl(AutoBuffer ab) {
      _gbCols = ab.getA8();
      int l = ab.get4();
      _agg = new AGG[l];
      for(int i=0;i<l;++i) _agg[i]=ab.get(AGG.class);
      int len = ab.get4();
      if( len == 0 ) return this;
      _g = new NonBlockingHashSet<>();
      for( int i=0;i<len;++i) _g.add(ab.get(G.class));
      return this;
    }
  }

  private static class GTask extends H2O.H2OCountedCompleter<GTask> {
    private final G _g;
    GTask(H2O.H2OCountedCompleter cc, G g) { super(cc); _g=g; }
    @Override protected void compute2() {
      _g.close();
      tryComplete();
    }
  }

  private static class ParallelPostGlobal extends H2O.H2OCountedCompleter<ParallelPostGlobal> {
    private final G[] _g;
    private final int _maxP=50*1000; // burn 50K at a time
    private final AtomicInteger _ctr;
    ParallelPostGlobal(G[] g) { _g=g; _ctr=new AtomicInteger(_maxP-1); }


    @Override protected void compute2(){
      addToPendingCount(_g.length-1);
      for( int i=0;i<Math.min(_g.length,_maxP);++i) frkTsk(i);
    }

    private void frkTsk(final int i) { new GTask(new Callback(), _g[i]).fork(); }

    private class Callback extends H2O.H2OCallback {
      public Callback(){super(ParallelPostGlobal.this);}
      @Override public void callback(H2O.H2OCountedCompleter cc) {
        int i = _ctr.incrementAndGet();
        if( i < _g.length )
          frkTsk(i);
      }
    }
  }

  private static class NBHSAD extends Iced {
    private transient NonBlockingHashSet _nd[];
    private int _n;
    NBHSAD(int n) { _nd = new NonBlockingHashSet[n]; _n=n; }
    @Override public AutoBuffer write_impl(AutoBuffer ab) {
      int len=_nd.length;
      ab.put4(len);
      for (NonBlockingHashSet a_nd : _nd) {
        if( a_nd==null ) {
          ab.put4(0);
          continue;
        }
        int s = a_nd.size();
        ab.put4(s);
        for (Object d : a_nd) ab.put8d((double)d);
      }
      return ab;
    }
    @Override public NBHSAD read_impl(AutoBuffer ab) {
      int len = ab.get4();
      _n=len;
      _nd=new NonBlockingHashSet[len];
      for(int i=0;i<len;++i) {
        _nd[i] = new NonBlockingHashSet<>();
        int s = ab.get4();
        if( s==0 ) continue;
        for(int j=0;j<s;++j) _nd[i].add(ab.get8d());
      }
      return this;
    }
  }

  private static class G extends Iced {

    public double _ds[];  // Array is final; contents change with the "fill"
    public int _hash;           // Hash is not final; changes with the "fill"
    public void fill(int row, Chunk chks[], long cols[]) {
      for( int c=0; c<cols.length; c++ ) // For all selection cols
        _ds[c] = chks[(int)cols[c]].atd(row); // Load into working array
      _hash = hash();
    }
    private int hash() {
      long h=0;                 // hash is sum of field bits
      for( double d : _ds ) h += Double.doubleToRawLongBits(d);
      // Doubles are lousy hashes; mix up the bits some
      h ^= (h>>>20) ^ (h>>>12);
      h ^= (h>>> 7) ^ (h>>> 4);
      return (int)((h^(h>>32))&0x7FFFFFFF);
    }
    @Override public boolean equals( Object o ) {
      return o instanceof G && Arrays.equals(_ds, ((G) o)._ds); }
    @Override public int hashCode() { return _hash; }
    @Override public String toString() { return Arrays.toString(_ds); }

    public long     _N;         // number of rows in the group, updated atomically
    public long[]   _ND;        // count of distincts, built from the NBHS<Double>
    public long[]   _NA;        // count of NAs for each aggregate, updated atomically
    public long[]   _f;         // first row, updated atomically
    public long[]   _l;         // last row, atomically updated
    public double[] _min;       // updated atomically
    public double[] _max;       // updated atomically
    public double[] _sum;       // sum, updated atomically
    public double[] _ss;        // sum of squares, updated atomically
    public double[] _avs;       // means, computed in the close
    public double[] _vars;      // vars, computed in the close
    public double[] _sdevs;     // sds,  computed in the close
//    private NBHSAD _nd;         // count distinct helper data structure
    private byte[] _NAMethod;

    // offset crud for unsafe
    private static final Unsafe U = UtilUnsafe.getUnsafe();
    private static final long _NOffset;

    // long[] offset and scale
    private static final int _8B = U.arrayBaseOffset(long[].class);
    private static final int _8S = U.arrayIndexScale(long[].class);
    // double[] offset and scale
    private static final int _dB = U.arrayBaseOffset(double[].class);
    private static final int _dS = U.arrayIndexScale(double[].class);

    // get the raw indices for the long[] and double[]
    private static long longRawIdx(int i)   { return _8B + _8S * i; }
    private static long doubleRawIdx(int i) { return _dB + _dS * i; }

    static {
      try {
        _NOffset   = U.objectFieldOffset(G.class.getDeclaredField("_N"));
      } catch(Exception e) { throw H2O.fail(); }
    }

    G(int row, Chunk[] cs, long[] cols,int aggs, byte[] naMethod) {
      _ds=new double[cols.length];
      this.fill(row, cs, cols);
      _NAMethod=naMethod;
//      _nd=new NBHSAD(aggs);
      _ND=new long[aggs];
      _NA=new long[aggs];
      _f =new long[aggs];
      _l =new long[aggs];
      _min=new double[aggs];
      _max=new double[aggs];
      _sum=new double[aggs];
      _ss =new double[aggs];
      _avs=new double[aggs];
      _vars=new double[aggs];
      _sdevs=new double[aggs];
      for( int i=0; i<_min.length; ++i) _min[i]=Double.POSITIVE_INFINITY;
      for( int i=0; i<_max.length; ++i) _max[i]=Double.NEGATIVE_INFINITY;
    }

    private void close() {
      for( int i=0;i<_NAMethod.length;++i ) {
        long n = _NAMethod[i]==AGG.T_RM?_N-_NA[i]:_N;
        _avs[i] = _sum[i]/n;
//        _ND[i] = _nd._nd[i]==null?0:_nd._nd[i].size(); _nd._nd[i]=null; // b free!
        _vars[i] = (_ss[i] - (_sum[i]*_sum[i])/n)/n;
        _sdevs[i]=Math.sqrt(_vars[i]);
      }
    }

    private static boolean CAS_N (G g, long o, long n          ) { return U.compareAndSwapLong(g,_NOffset,o,n); }
    private static boolean CAS_NA(G g, long off, long o, long n) { return U.compareAndSwapLong(g._NA,off,o,n);  }
    private static boolean CAS_f (G g, long off, long o, long n) { return U.compareAndSwapLong(g._f,off,o,n);   }
    private static boolean CAS_l (G g, long off, long o, long n) { return U.compareAndSwapLong(g._l,off,o,n);   }

    // doubles are toRawLongBits'ized, and passed as longs
    private static boolean CAS_min(G g, long off, long o, long n) { return U.compareAndSwapLong(g._min,off,o,n);}
    private static boolean CAS_max(G g, long off, long o, long n) { return U.compareAndSwapLong(g._max,off,o,n);}
    private static boolean CAS_sum(G g, long off, long o, long n) { return U.compareAndSwapLong(g._sum,off,o,n);}
    private static boolean CAS_ss (G g, long off, long o, long n) { return U.compareAndSwapLong(g._ss ,off,o,n);}
  }

  static class AGG extends AST {
    // (AGG #N "agg" #col "na"  "agg" #col "na"   => string num string   string num string
    private AGG[] _aggs;
    protected AGG parse_impl(Exec E) {
      int n = (int)((ASTNum)(E.parse()))._d; E.skipWS();
      _aggs=new AGG[n];
      for( int i=0;i<n;++i) {
        String type = E.parseString(E.peekPlus()); E.skipWS();
        int     col = (int)((ASTNum)E.parse()).dbl(); E.skipWS();
        String   na = E.parseString(E.peekPlus()); E.skipWS();
        String name = E.parseString(E.peekPlus()); E.skipWS();
        _aggs[i]=new AGG(type,col,na,name);
      }
      return this;
    }

    // Aggregate types
    private static final byte T_N  = 0;
    private static final byte T_ND = 1;
    private static final byte T_F  = 2;
    private static final byte T_L  = 3;
    private static final byte T_MIN= 4;
    private static final byte T_MAX= 5;
    private static final byte T_AVG= 6;
    private static final byte T_SD = 7;
    private static final byte T_VAR= 8;
    private static final byte T_SUM= 9;
    private static final byte T_SS = 10;

    // How to handle NAs
    private static final byte T_ALL = 0;
    private static final byte T_IG  = 1;
    private static final byte T_RM  = 2;

    private static transient HashMap<String,Byte> TM = new HashMap<>();
    static{
      // aggregates
      TM.put("count",       (byte)0);
      TM.put("count_unique",(byte)1);
      TM.put("first",       (byte)2);
      TM.put("last",        (byte)3);
      TM.put("min",         (byte)4);
      TM.put("max",         (byte)5);
      TM.put("mean",        (byte)6);
      TM.put("avg",         (byte)6);
      TM.put("sd",          (byte)7);
      TM.put("stdev",       (byte)7);
      TM.put("var",         (byte)8);
      TM.put("sum",         (byte)9);
      TM.put("ss",          (byte)10);
      // na handling
      TM.put("ignore"      ,(byte)0);
      TM.put("rm"          ,(byte)1);
      TM.put("all"         ,(byte)2);
    }

    private final byte _type;
    private final int _c;
    private final String _name;
    private final byte _na_handle;
    AGG() {_type=0;_c=-1;_name=null;_na_handle=0;}
    AGG(String s, int c, String na, String name) {
      _type=TM.get(s.toLowerCase());
      _c=c;
      _name=(name==null || name.equals(""))?s+"_C"+(c+1):name;
      if( !TM.keySet().contains(na) ) {
        Log.info("Unknown NA handle type given: `" + na + "`. Switching to \"ignore\" method.");
        _na_handle=0;
      } else _na_handle = TM.get(na);
    }

    private static String[] names(AGG[] _agg) {
      String[] names = new String[_agg.length];
      for(int i=0;i<names.length;++i)
        names[i] = _agg[i]._name;
      return names;
    }

    private static byte[] naMethods(AGG[] agg) {
      byte[] methods = new byte[agg.length];
      for(int i=0;i<agg.length;++i)
        methods[i]=agg[i]._na_handle;
      return methods;
    }

    private boolean isIgnore() { return _na_handle == 0; }
    private boolean isRemove() { return _na_handle == 1; }
    private boolean isAll()    { return _na_handle == 2; }

    // satisfy the extends
    @Override void exec(Env e) { throw H2O.fail();}
    @Override String value() { return "agg"; }
    @Override int type() { return 0; }
  }
}

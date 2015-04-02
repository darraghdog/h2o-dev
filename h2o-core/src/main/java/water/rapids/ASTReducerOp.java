package water.rapids;

import java.util.ArrayList;
import water.*;
import water.fvec.*;

// Variable length; instances will be created of required length
public abstract class ASTReducerOp extends ASTOp {
  protected static double _init;
  protected static boolean _narm;        // na.rm in R
  protected static int _argcnt;
  protected ASTReducerOp( double init) {
    super(new String[]{"","dblary","...", "na.rm"});
    _init = init;
  }

  protected ASTReducerOp parse_impl(Exec E) {
    ArrayList<AST> dblarys = new ArrayList<>();
    AST ary = E.parse();
    if (ary instanceof ASTId) ary = Env.staticLookup((ASTId)ary);
    dblarys.add(ary);
    AST a;
    E.skipWS();
    while (true) {
      a = E.skipWS().parse();
      if (a instanceof ASTId) {
        AST ast = E._env.lookup((ASTId)a);
        if (ast instanceof ASTFrame) {dblarys.add(a); continue; } else break;
      }
      if (a instanceof ASTNum || a instanceof ASTFrame || a instanceof ASTSlice || a instanceof ASTBinOp || a instanceof ASTUniOp || a instanceof ASTReducerOp)
        dblarys.add(a);
      else break;
    }
    // Get the na.rm last
    try {
      a = E._env.lookup((ASTId) a);
    } catch (ClassCastException e) {
      throw new IllegalArgumentException("Expected the na.rm value to be one of $TRUE, $FALSE, $T, $F");
    }
    _narm = ((ASTNum)a).dbl() == 1;
    ASTReducerOp res = (ASTReducerOp) clone();
    AST[] arys = new AST[_argcnt = dblarys.size()];
    for (int i = 0; i < dblarys.size(); i++) arys[i] = dblarys.get(i);
    res._asts = arys;
    return res;
  }

  @Override double[] map(Env env, double[] in, double[] out, AST[] args) {
    double s = _init;
    for (double v : in) if (!_narm || !Double.isNaN(v)) s = op(s,v);
    if (out == null || out.length < 1) out = new double[1];
    out[0] = s;
    return out;
  }
  protected abstract double op( double d0, double d1 );
  @Override protected void apply(Env env) {
    double sum=_init;
    int argcnt = _argcnt;
    for( int i=0; i<argcnt; i++ )
      if( env.isNum() ) sum = op(sum,env.popDbl());
      else {
        Frame fr = env.popAry(); // pop w/o lowering refcnts ... clean it up later
        for(Vec v : fr.vecs()) if (v.isEnum() || v.isUUID() || v.isString()) throw new IllegalArgumentException("`"+opStr()+"`" + " only defined on a data frame with all numeric variables");
        sum = op(sum,_narm?new NaRmRedOp(this).doAll(fr)._d:new RedOp(this).doAll(fr)._d);
      }
    env.push(new ValNum(sum));
  }

  @Override void exec(Env e, AST arg1, AST[] args) {
    if (args == null) {
      _init = 0;
      _narm = true;
      _argcnt = 1;
    }
    arg1.exec(e);
    e._global._frames.put(Key.make().toString(), e.peekAry());
    apply(e);
  }

  private static class RedOp extends MRTask<RedOp> {
    final ASTReducerOp _bin;
    RedOp( ASTReducerOp bin ) { _bin = bin; _d = ASTReducerOp._init; }
    double _d;
    @Override public void map( Chunk chks[] ) {
      int rows = chks[0]._len;
      for (Chunk C : chks) {
        assert C.vec().isNumeric();
        double sum = _d;
        for (int r = 0; r < rows; r++)
          sum = _bin.op(sum, C.atd(r));
        _d = sum;
        if (Double.isNaN(sum)) break;
      }
    }
    @Override public void reduce( RedOp s ) { _d = _bin.op(_d,s._d); }
  }

  private static class NaRmRedOp extends MRTask<NaRmRedOp> {
    final ASTReducerOp _bin;
    NaRmRedOp( ASTReducerOp bin ) { _bin = bin; _d = ASTReducerOp._init; }
    double _d;
    @Override public void map( Chunk chks[] ) {
      int rows = chks[0]._len;
      for (Chunk C : chks) {
        assert C.vec().isNumeric();
        double sum = _d;
        for (int r = 0; r < rows; r++) {
          double d = C.atd(r);
          if (!Double.isNaN(d))
            sum = _bin.op(sum, d);
        }
        _d = sum;
        if (Double.isNaN(sum)) break;
      }
    }
    @Override public void reduce( NaRmRedOp s ) { _d = _bin.op(_d,s._d); }
  }
}

class ASTSum extends ASTReducerOp { 
  ASTSum() {super(0);} 
  @Override protected String opStr(){ return "sum";} 
  @Override protected ASTOp make() {return new ASTSum();} 
  @Override protected double op(double d0, double d1) { return d0+d1;}
  @Override protected void apply(Env env) {
    double sum=_init;
    int argcnt = _argcnt;
    for( int i=0; i<argcnt; i++ )
      if( env.isNum() ) sum = op(sum,env.popDbl());
      else {
        Frame fr = env.popAry(); // pop w/o lowering refcnts ... clean it up later
        for(Vec v : fr.vecs()) if (v.isEnum() || v.isUUID() || v.isString()) throw new IllegalArgumentException("`"+opStr()+"`" + " only defined on a data frame with all numeric variables");
        sum += new RedSum(_narm).doAll(fr)._d;
      }
    env.push(new ValNum(sum));
  }

  private static class RedSum extends MRTask<RedSum> {
    final boolean _narm;
    double _d;
    RedSum( boolean narm ) { _narm = narm; }
    @Override public void map( Chunk chks[] ) {
      int rows = chks[0]._len;
      for (Chunk C : chks) {
        assert C.vec().isNumeric();
        double sum=_d;
        if( _narm ) for (int r = 0; r < rows; r++) { double d = C.atd(r); if( !Double.isNaN(d) ) sum += d; }
        else        for (int r = 0; r < rows; r++) { double d = C.atd(r);                        sum += d; }
        _d = sum;
        if( Double.isNaN(sum) ) break;
      }
    }
    @Override public void reduce( RedSum s ) { _d += s._d; }
  }
}


class ASTMin extends ASTReducerOp {
  ASTMin( ) { super( Double.POSITIVE_INFINITY); }
  @Override protected String opStr(){ return "min";}
  @Override protected ASTOp make() {return new ASTMin();}
  @Override protected double op(double d0, double d1) { return Math.min(d0, d1); }
  protected ASTMin parse_impl(Exec E) { return (ASTMin)super.parse_impl(E); }
  @Override protected void apply(Env env) {
    double min = Double.POSITIVE_INFINITY;
    int argcnt = env.sp();
    for( int i=0; i<argcnt; i++ )
      if( env.isNum() ) min = Math.min(min, env.popDbl());
      else {
        Frame fr = env.popAry();
        for(Vec v : fr.vecs()) if (v.isEnum() || v.isUUID() || v.isString()) throw new IllegalArgumentException("`"+opStr()+"`" + " only defined on a data frame with all numeric variables");
        for (Vec v : fr.vecs())
          if (v.naCnt() > 0 && !_narm) { min = Double.NaN; break; }
          else min = Math.min(min, v.min());
      }
    env.push(new ValNum(min));
  }
}

class ASTMax extends ASTReducerOp {
  ASTMax( ) { super( Double.NEGATIVE_INFINITY); }
  @Override protected String opStr(){ return "max";}
  @Override protected ASTOp make() {return new ASTMax();}
  @Override protected double op(double d0, double d1) { return Math.max(d0,d1); }
  protected ASTMax parse_impl(Exec E) { return (ASTMax)super.parse_impl(E); }
  @Override protected void apply(Env env) {
    double max = Double.NEGATIVE_INFINITY;
    int argcnt = env.sp();
    for( int i=0; i<argcnt; i++ )
      if( env.isNum() ) max = Math.max(max, env.popDbl());
      else {
        Frame fr = env.popAry();
        for(Vec v : fr.vecs()) if (v.isEnum() || v.isUUID() || v.isString()) throw new IllegalArgumentException("`"+opStr()+"`" + " only defined on a data frame with all numeric variables");
        for (Vec v : fr.vecs())
          if (v.naCnt() > 0 && !_narm) { max = Double.NaN; break; }
          else max = Math.max(max, v.max());
      }
    env.push(new ValNum(max));
  }
}


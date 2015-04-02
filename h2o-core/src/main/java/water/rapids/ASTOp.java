package water.rapids;

import hex.DMatrix;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import org.joda.time.DateTime;
import org.joda.time.format.DateTimeFormatter;
import water.*;
import water.fvec.*;
import water.nbhm.NonBlockingHashMap;
import water.parser.ParseTime;
import water.parser.ValueString;
import water.util.MathUtils;

/**
 * Parse a generic R string and build an AST, in the context of an H2O Cloud
 */

// --------------------------------------------------------------------------
public abstract class ASTOp extends AST {
  // Tables of operators by arity
  static final public HashMap<String,ASTOp> UNI_INFIX_OPS  = new HashMap<>();
  static final public HashMap<String,ASTOp> BIN_INFIX_OPS  = new HashMap<>();
  static final public HashMap<String,ASTOp> PREFIX_OPS     = new HashMap<>();
  static final public HashMap<String,ASTOp> UDF_OPS        = new HashMap<>();
  static final public HashMap<String, AST>  SYMBOLS        = new HashMap<>();
  // Too avoid a cyclic class-loading dependency, these are init'd before subclasses.
  static final String VARS1[] = new String[]{ "", "x"};
  static final String VARS2[] = new String[]{ "", "x","y"};
  static {
    // All of the special chars (see Exec.java)
    SYMBOLS.put(",", new ASTStatement());
    SYMBOLS.put("=", new ASTAssign());
    SYMBOLS.put("'", new ASTString('\'', ""));
    SYMBOLS.put("\"",new ASTString('\"', ""));
    SYMBOLS.put("%", new ASTId('%', ""));
    SYMBOLS.put("!", new ASTId('!', ""));
    SYMBOLS.put("#", new ASTNum(0));
    SYMBOLS.put("g", new ASTGT());
    SYMBOLS.put("G", new ASTGE());
    SYMBOLS.put("l", new ASTLT());
    SYMBOLS.put("L", new ASTLE());
    SYMBOLS.put("N", new ASTNE());
    SYMBOLS.put("n", new ASTEQ());
    SYMBOLS.put("[", new ASTSlice());
    SYMBOLS.put("{", new ASTSeries(null, null));
    SYMBOLS.put(":", new ASTSpan(new ASTNum(0),new ASTNum(0)));
    SYMBOLS.put("not", new ASTNot());
    SYMBOLS.put("_", new ASTNot());
    SYMBOLS.put("if", new ASTIf());
    SYMBOLS.put("else", new ASTElse());
    SYMBOLS.put("for", new ASTFor());
    SYMBOLS.put("while", new ASTWhile());
    SYMBOLS.put("return", new ASTReturn());
    SYMBOLS.put("del", new ASTDelete());
    SYMBOLS.put("x", new ASTMMult());
    SYMBOLS.put("t", new ASTTranspose());
    SYMBOLS.put("agg",new ASTGroupBy.AGG());

    //TODO: Have `R==` type methods (also `py==`, `js==`, etc.)

    // Unary infix ops
    putUniInfix(new ASTNot());
    // Binary infix ops
    putBinInfix(new ASTPlus());
    putBinInfix(new ASTSub());
    putBinInfix(new ASTMul());
    putBinInfix(new ASTMMult());
    putBinInfix(new ASTDiv());
    putBinInfix(new ASTPow());
    putBinInfix(new ASTPow2());
    putBinInfix(new ASTMod());
    putBinInfix(new ASTAND());
    putBinInfix(new ASTOR());
    putBinInfix(new ASTLT());
    putBinInfix(new ASTLE());
    putBinInfix(new ASTGT());
    putBinInfix(new ASTGE());
    putBinInfix(new ASTEQ());
    putBinInfix(new ASTNE());
    putBinInfix(new ASTLA());
    putBinInfix(new ASTLO());

    // Unary prefix ops
    putPrefix(new ASTIsNA());
    putPrefix(new ASTNrow());
    putPrefix(new ASTNcol());
    putPrefix(new ASTLength());
    putPrefix(new ASTAbs ());
    putPrefix(new ASTSgn ());
    putPrefix(new ASTSqrt());
    putPrefix(new ASTCeil());
    putPrefix(new ASTFlr ());
    putPrefix(new ASTLog ());
    putPrefix(new ASTLog10 ());
    putPrefix(new ASTLog2 ());
    putPrefix(new ASTLog1p ());
    putPrefix(new ASTExp ());
    putPrefix(new ASTExpm1 ());
    putPrefix(new ASTGamma());
    putPrefix(new ASTLGamma());
    putPrefix(new ASTDiGamma());
    putPrefix(new ASTTriGamma());
    putPrefix(new ASTScale());
    putPrefix(new ASTCharacter());
    putPrefix(new ASTFactor());
    putPrefix(new ASTAsNumeric());
    putPrefix(new ASTIsFactor());
    putPrefix(new ASTAnyFactor());              // For Runit testing
    putPrefix(new ASTCanBeCoercedToLogical());
    putPrefix(new ASTAnyNA());
    putPrefix(new ASTRound());
    putPrefix(new ASTSignif());
    putPrefix(new ASTTrun());

    putPrefix(new ASTTranspose());

    // Trigonometric functions
    putPrefix(new ASTCos  ());
    putPrefix(new ASTSin  ());
    putPrefix(new ASTTan  ());
    putPrefix(new ASTACos ());
    putPrefix(new ASTASin ());
    putPrefix(new ASTATan ());
    putPrefix(new ASTCosh ());
    putPrefix(new ASTSinh ());
    putPrefix(new ASTTanh ());
    putPrefix(new ASTACosh());
    putPrefix(new ASTASinh());
    putPrefix(new ASTATanh());
    putPrefix(new ASTCosPi());
    putPrefix(new ASTSinPi());
    putPrefix(new ASTTanPi());


    // More generic reducers
    putPrefix(new ASTMin ());
    putPrefix(new ASTMax ());
    putPrefix(new ASTSum ());
    putPrefix(new ASTSdev());
    putPrefix(new ASTVar ());
    putPrefix(new ASTMean());

    // Misc
    putPrefix(new ASTSetLevel());
    putPrefix(new ASTMatch ());
    putPrefix(new ASTRename());
    putPrefix(new ASTSeq   ());
    putPrefix(new ASTSeqLen());
    putPrefix(new ASTRepLen());
    putPrefix(new ASTCbind ());
    putPrefix(new ASTRbind ());
    putPrefix(new ASTTable ());
//    putPrefix(new ASTReduce());
    putPrefix(new ASTIfElse());
    putPrefix(new ASTApply ());
    putPrefix(new ASTSApply());
    putPrefix(new ASTddply());
    putPrefix(new ASTMerge ());
    putPrefix(new ASTGroupBy());
//    putPrefix(new ASTUnique());
    putPrefix(new ASTXorSum());
    putPrefix(new ASTRunif ());
    putPrefix(new ASTCut   ());
    putPrefix(new ASTLs    ());
    putPrefix(new ASTSetColNames());

    // Date
    putPrefix(new ASTasDate());

//Classes that may not come back:

//    putPrefix(new ASTfindInterval());
//    putPrefix(new ASTPrint ());
    putPrefix(new ASTCat   ());
//Time extractions, to and from msec since the Unix Epoch
    putPrefix(new ASTYear  ());
    putPrefix(new ASTMonth ());
    putPrefix(new ASTDay   ());
    putPrefix(new ASTDayOfWeek());
    putPrefix(new ASTHour  ());
    putPrefix(new ASTMinute());
    putPrefix(new ASTSecond());
    putPrefix(new ASTMillis());
    putPrefix(new ASTMktime());
    putPrefix(new ASTFoldCombine());
    putPrefix(new COp());
    putPrefix(new ROp());
    putPrefix(new O());

//    // Time series operations
//    putPrefix(new ASTDiff  ());
//    putPrefix(new ASTIsTRUE());
//    putPrefix(new ASTMTrans());
//    putBinInfix(new ASTMMult());
  }
  static public void registerAST(ASTOp op) { putPrefix(op); }

  static private void putUniInfix(ASTOp ast) { UNI_INFIX_OPS.put(ast.opStr(),ast); }
  static private void putBinInfix(ASTOp ast) { BIN_INFIX_OPS.put(ast.opStr(),ast); SYMBOLS.put(ast.opStr(), ast); }
  static private void putPrefix  (ASTOp ast) { PREFIX_OPS.put(ast.opStr(),ast);    SYMBOLS.put(ast.opStr(), ast); }
  static         void putUDF     (ASTOp ast, String fn) { UDF_OPS.put(fn, ast);    SYMBOLS.put(fn, ast);}
  static public boolean isUDF(String id) { return UDF_OPS.containsKey(id); }

  // All fields are final, because functions are immutable
  final String _vars[]; // Variable names
  ASTOp( String vars[]) { _vars = vars; }

  abstract protected String opStr();
  abstract protected  ASTOp  make();
  // Standard column-wise function application
  abstract protected void apply(Env e);
  // Special row-wise 'apply'
  double[] map(Env env, double[] in, double[] out, AST[] args) { throw H2O.unimpl(); }
  @Override void exec(Env e) { throw H2O.fail(); }
  // special exec for apply calls
  void exec(Env e, AST arg1, AST[] args) { throw H2O.unimpl("No exec method for `" + this.opStr() + "` during `apply` call"); }
  @Override int type() { throw H2O.fail(); }
  @Override String value() { throw H2O.fail(); }

  public static ASTOp get(String op) {
    if (BIN_INFIX_OPS.containsKey(op)) return BIN_INFIX_OPS.get(op);
    if (UNI_INFIX_OPS.containsKey(op)) return UNI_INFIX_OPS.get(op);
    if (isUDF(op)) return UDF_OPS.get(op);
    if (PREFIX_OPS.containsKey(op)) return PREFIX_OPS.get(op);
    throw H2O.fail("Unimplemented: Could not find the operation or function "+op);
  }
}

abstract class ASTUniOrBinOp extends ASTOp {
  ASTUniOrBinOp(String[] vars) { super(vars); }
  double op( double d ) { throw H2O.fail(); }
  double op(double d0, double d1) { throw H2O.fail(); }
  String op( String s0, double d1 ) { throw H2O.fail(); }
  String op( double d0, String s1 ) { throw H2O.fail(); }
  String op( String s0, String s1 ) { throw H2O.fail(); }
}

abstract class ASTUniOp extends ASTUniOrBinOp {
  ASTUniOp() { super(VARS1); }
  protected ASTUniOp( String[] vars) { super(vars); }
  protected ASTUniOp parse_impl(Exec E) {
    if (!E.hasNext()) throw new IllegalArgumentException("End of input unexpected. Badly formed AST.");
    AST arg = E.parse();
    if (arg instanceof ASTId) arg = Env.staticLookup((ASTId)arg);
    ASTUniOp res = (ASTUniOp) clone();
    res._asts = new AST[]{arg};
    return res;
  }

  @Override void exec(Env e, AST arg1, AST[] args) {
    if (args != null) throw new IllegalArgumentException("Too many arguments passed to `"+opStr()+"`");
    arg1.exec(e);
    if (e.isAry()) e._global._frames.put(Key.make().toString(), e.peekAry());
    apply(e);
  }

  @Override protected void apply(Env env) {
    // Expect we can broadcast across all functions as needed.
    if( env.isNum() ) { env.push(new ValNum(op(env.popDbl()))); return; }
//    if( env.isStr() ) { env.push(new ASTString(op(env.popStr()))); return; }
    Frame fr = env.popAry();
    final ASTUniOp uni = this;  // Final 'this' so can use in closure
    Frame fr2 = new MRTask() {
      @Override public void map( Chunk[] chks, NewChunk[] nchks ) {
        for( int i=0; i<nchks.length; i++ ) {
          NewChunk n =nchks[i];
          Chunk c = chks[i];
          int rlen = c._len;
          if (c.vec().isEnum() || c.vec().isUUID() || c.vec().isString()) {
            for (int r = 0; r <rlen;r++) n.addNum(Double.NaN);
          } else {
            for( int r=0; r<rlen; r++ )
              n.addNum(uni.op(c.atd(r)));
          }
        }
      }
    }.doAll(fr.numCols(),fr).outputFrame(fr._names, null);
    env.pushAry(fr2);
  }
}

class ASTasDate extends ASTOp {
  protected static String _format;
  ASTasDate() { super(new String[]{"as.Date", "x", "format"}); }
  @Override protected String opStr() { return "as.Date"; }
  @Override protected ASTOp make() {return new ASTasDate();}
  @Override
  protected ASTasDate parse_impl(Exec E) {
    AST ast = E.parse();
    if (ast instanceof ASTId) ast = Env.staticLookup((ASTId)ast);
    try {
      _format = ((ASTString)E.skipWS().parse())._s;
    } catch (ClassCastException e) {
      throw new IllegalArgumentException("`format` must be a string.");
    }
    ASTasDate res = (ASTasDate) clone();
    res._asts = new AST[]{ast};
    return res;
  }
  @Override protected void apply(Env env) {
    final String format = _format;
    if (format.isEmpty()) throw new IllegalArgumentException("as.Date requires a non-empty format string");
    // check the format string more?

    Frame fr = env.popAry();

    if( fr.vecs().length != 1 || !fr.vecs()[0].isEnum() )
      throw new IllegalArgumentException("as.Date requires a single column of factors");

    Frame fr2 = new MRTask() {
      @Override public void map( Chunk chks[], NewChunk nchks[] ) {
        //done on each node in lieu of rewriting DateTimeFormatter as Iced
        DateTimeFormatter dtf = ParseTime.forStrptimePattern(format).withZone(ParseTime.getTimezone());
        for( int i=0; i<nchks.length; i++ ) {
          NewChunk n =nchks[i];
          Chunk c = chks[i];
          int rlen = c._len;
          for( int r=0; r<rlen; r++ ) {
            if (!c.isNA(r)) {
              String date = c.vec().domain()[(int)c.atd(r)];
              n.addNum(DateTime.parse(date, dtf).getMillis(), 0);
            } else n.addNA();
          }
        }
      }
    }.doAll(fr.numCols(),fr).outputFrame(fr._names, null);
    env.pushAry(fr2);
  }
}

//
//// Finite backward difference for user-specified lag
//// http://en.wikipedia.org/wiki/Finite_difference
//class ASTDiff extends ASTOp {
//  ASTDiff() { super(new String[]{"diff", "x", "lag", "differences"},
//          new Type[]{Type.ARY, Type.ARY, Type.DBL, Type.DBL},
//          OPF_PREFIX,
//          OPP_PREFIX,
//          OPA_RIGHT); }
//  @Override protected String opStr() { return "diff"; }
//  @Override protected ASTOp make() {return new ASTDiff();}
//  @Override protected void apply(Env env, int argcnt, ASTApply apply) {
//    final int diffs = (int)env.popDbl();
//    if(diffs < 0) throw new IllegalArgumentException("differences must be an integer >= 1");
//    final int lag = (int)env.popDbl();
//    if(lag < 0) throw new IllegalArgumentException("lag must be an integer >= 1");
//
//    Frame fr = env.popAry();
//    String skey = env.key();
//    if(fr.vecs().length != 1 || fr.vecs()[0].isEnum())
//      throw new IllegalArgumentException("diff takes a single numeric column vector");
//
//    Frame fr2 = new MRTask() {
//      @Override public void map(Chunk chk, NewChunk nchk) {
//        int rstart = (int)(diffs*lag - chk._start);
//        if(rstart > chk._len) return;
//        rstart = Math.max(0, rstart);
//
//        // Formula: \Delta_h^n x_t = \sum_{i=0}^n (-1)^i*\binom{n}{k}*x_{t-i*h}
//        for(int r = rstart; r < chk._len; r++) {
//          double x = chk.atd(r);
//          long row = chk._start + r;
//
//          for(int i = 1; i <= diffs; i++) {
//            double x_lag = chk.at_slow(row - i*lag);
//            double coef = ArithmeticUtils.binomialCoefficient(diffs, i);
//            x += (i % 2 == 0) ? coef*x_lag : -coef*x_lag;
//          }
//          nchk.addNum(x);
//        }
//      }
//    }.doAll(1,fr).outputFrame(fr.names(), fr.domains());
//    env.subRef(fr, skey);
//    env.pop();
//    env.push(fr2);
//  }
//}


/**
 *  ASTBinOp: E x E -&gt; E
 *
 *  This covers the class of operations that produce an array, scalar, or string from the cartesian product
 *  of the set E = {x | x is a string, scalar, or array}.
 */
abstract class ASTBinOp extends ASTUniOrBinOp {

  ASTBinOp() { super(VARS2); } // binary ops are infix ops

  protected ASTBinOp parse_impl(Exec E) {
    AST l = E.parse();
    if (l instanceof ASTId) l = Env.staticLookup((ASTId)l);
    AST r = E.parse();
    if (r instanceof ASTId) r = Env.staticLookup((ASTId)r);
    ASTBinOp res = (ASTBinOp) clone();
    res._asts = new AST[]{l,r};
    return res;
  }

  @Override protected void apply(Env env) {
    // Expect we can broadcast across all functions as needed.
    Frame fr0 = null, fr1 = null;
    double d0=0, d1=0;
    String s0=null, s1=null;

    // Must pop ONLY twice off the stack
    int left_type = env.peekType();
    Object left = env.peek();
    int right_type = env.peekTypeAt(-1);
    Object right = env.peekAt(-1);

    // Cast the LHS of the op
    switch(left_type) {
      case Env.NUM: d0  = ((ValNum)left)._d; break;
      case Env.ARY: fr0 = ((ValFrame)left)._fr; break;
      case Env.STR: s0  = ((ValStr)left)._s; break;
      default: throw H2O.fail("Got unusable type: "+ left_type +" in binary operator "+ opStr());
    }

    // Cast the RHS of the op
    switch(right_type) {
      case Env.NUM: d1  = ((ValNum)right)._d; break;
      case Env.ARY: fr1 = ((ValFrame)right)._fr; break;
      case Env.STR: s1  = ((ValStr)right)._s; break;
      default: throw H2O.fail("Got unusable type: "+ right_type +" in binary operator "+ opStr());
    }

    // If both are doubles on the stack
    if( (fr0==null && fr1==null) && (s0==null && s1==null) ) { env.pop(); env.pop(); env.push(new ValNum(op(d0, d1))); return; }

    // One or both of the items on top of stack are Strings and neither are frames
    if( fr0==null && fr1==null) {
      env.pop(); env.pop();
      // s0 == null -> op(d0, s1)
      if (s0 == null) {
        // cast result of op if doing comparison, else combine the Strings if defined for op
        if (opStr().equals("==") || opStr().equals("!=")) env.push(new ValNum(Double.valueOf(op(d0,s1))));
        else env.push(new ValStr(op(d0,s1)));
      }
      // s1 == null -> op(s0, d1)
      else if (s1 == null) {
        // cast result of op if doing comparison, else combine the Strings if defined for op
        if (opStr().equals("==") || opStr().equals("!=")) env.push(new ValNum(Double.valueOf(op(s0,d1))));
        else env.push(new ValStr(op(s0,d1)));
      // s0 != null, s1 != null
      } else env.push(new ValStr(op(s0,s1)));
      return;
    }

    final boolean lf = fr0 != null;
    final boolean rf = fr1 != null;
    final double df0 = d0, df1 = d1;
    final String sf0 = s0, sf1 = s1;
    Frame fr;           // Do-All frame
    int ncols;          // Result column count
    if( fr0 !=null ) {  // Left?
      ncols = fr0.numCols();
      if( fr1 != null ) {
        if( fr0.numCols() != fr1.numCols() ||
            fr0.numRows() != fr1.numRows() )
          throw new IllegalArgumentException("Arrays must be same size: LHS FRAME NUM ROWS/COLS: "+fr0.numRows()+"/"+fr0.numCols() +" vs RHS FRAME NUM ROWS/COLS: "+fr1.numRows()+"/"+fr1.numCols());
        fr = new Frame(fr0).add(fr1);
      } else {
        fr = new Frame(fr0);
      }
    } else {
      ncols = fr1.numCols();
      fr = new Frame(fr1);
    }
    final ASTBinOp bin = this;  // Final 'this' so can use in closure

    // Run an arbitrary binary op on one or two frames & scalars
    Frame fr2 = new MRTask() {
      @Override public void map( Chunk chks[], NewChunk nchks[] ) {
        for( int i=0; i<nchks.length; i++ ) {
          NewChunk n =nchks[i];
          int rlen = chks[0]._len;
          Chunk c0 = chks[i];
          if( (!c0.vec().isEnum() &&
                  !(lf && rf && chks[i+nchks.length].vec().isEnum())) ||
                  bin instanceof ASTEQ ||
                  bin instanceof ASTNE ) {

            // Loop over rows
            for( int ro=0; ro<rlen; ro++ ) {
              double lv=0; double rv=0; String l=null; String r=null;

              // Initialize the lhs value
              if (lf) {
                if(chks[i].vec().isUUID() || (chks[i].isNA(ro) && !bin.opStr().equals("|"))) { n.addNum(Double.NaN); continue; }
                if (chks[i].vec().isEnum()) l = chks[i].vec().domain()[(int)chks[i].atd(ro)];
                else lv = chks[i].atd(ro);
              } else if (sf0 == null) {
                if (Double.isNaN(df0) && !bin.opStr().equals("|")) { n.addNum(Double.NaN); continue; }
                lv = df0; l = null;
              } else {
                l = sf0;
              }

              // Initialize the rhs value
              if (rf) {
                if(chks[i+(lf ? nchks.length:0)].vec().isUUID() || chks[i].isNA(ro) && !bin.opStr().equals("|")) { n.addNum(Double.NaN); continue; }
                if (chks[i].vec().isEnum()) r = chks[i].vec().domain()[(int)chks[i].atd(ro)];
                else rv = chks[i+(lf ? nchks.length:0)].atd(ro);
              } else if (sf1 == null) {
                if (Double.isNaN(df1) && !bin.opStr().equals("|")) { n.addNum(Double.NaN); continue; }
                rv = df1; r= null;
              } else {
                r = sf1;
              }

              // Append the value to the chunk after applying op(lhs,rhs)
              if (l == null && r == null)
                n.addNum(bin.op(lv, rv));
              else if (l == null) n.addNum(Double.valueOf(bin.op(lv,r)));
              else if (r == null) n.addNum(Double.valueOf(bin.op(l,rv)));
              else n.addNum(Double.valueOf(bin.op(l,r)));
            }
          } else {
            for( int r=0; r<rlen; r++ )  n.addNA();
          }
        }
      }
    }.doAll(ncols,fr).outputFrame(null, (lf ? fr0 : fr1)._names,null);
    env.poppush(2, new ValFrame(fr2));
  }
  @Override public String toString() { return "("+opStr()+" "+Arrays.toString(_asts)+")"; }
}

class ASTPlus extends ASTBinOp { public ASTPlus() { super(); } @Override protected String opStr(){ return "+";} @Override protected ASTOp make() {return new ASTPlus();}
  @Override double op(double d0, double d1) { return d0+d1;}
  @Override String op(String s0, double d1) {throw new IllegalArgumentException("Cannot add Strings.");}
  @Override String op(double d0, String s1) {throw new IllegalArgumentException("Cannot add Strings.");}
  @Override String op(String s0, String s1) {throw new IllegalArgumentException("Cannot add Strings.");}
}
class ASTSub extends ASTBinOp { public ASTSub() { super(); } @Override protected String opStr(){ return "-";} @Override protected ASTOp make() {return new ASTSub ();}
  @Override double op(double d0, double d1) { return d0-d1;}
  @Override String op(String s0, double d1) {throw new IllegalArgumentException("Cannot subtract Strings.");}
  @Override String op(double d0, String s1) {throw new IllegalArgumentException("Cannot subtract Strings.");}
  @Override String op(String s0, String s1) {throw new IllegalArgumentException("Cannot subtract Strings.");}
}
class ASTMul extends ASTBinOp { public ASTMul() { super(); } @Override protected String opStr(){ return "*";} @Override protected ASTOp make() {return new ASTMul ();}
  @Override double op(double d0, double d1) { return d0*d1;}
  @Override String op(String s0, double d1) {throw new IllegalArgumentException("Cannot multiply Strings.");}
  @Override String op(double d0, String s1) {throw new IllegalArgumentException("Cannot multiply Strings.");}
  @Override String op(String s0, String s1) {throw new IllegalArgumentException("Cannot multiply Strings.");}
}
class ASTDiv extends ASTBinOp { public ASTDiv() { super(); } @Override protected String opStr(){ return "/";} @Override protected ASTOp make() {return new ASTDiv ();}
  @Override double op(double d0, double d1) { return d0/d1;}
  @Override String op(String s0, double d1) {throw new IllegalArgumentException("Cannot divide Strings.");}
  @Override String op(double d0, String s1) {throw new IllegalArgumentException("Cannot divide Strings.");}
  @Override String op(String s0, String s1) {throw new IllegalArgumentException("Cannot divide Strings.");}
}
class ASTPow extends ASTBinOp { public ASTPow() { super(); } @Override protected String opStr(){ return "^"  ;} @Override protected ASTOp make() {return new ASTPow ();}
  @Override double op(double d0, double d1) { return Math.pow(d0,d1);}
  @Override String op(String s0, double d1) {throw new IllegalArgumentException("Cannot exponentiate Strings.");}
  @Override String op(double d0, String s1) {throw new IllegalArgumentException("Cannot exponentiate Strings.");}
  @Override String op(String s0, String s1) {throw new IllegalArgumentException("Cannot exponentiate Strings.");}
}
class ASTPow2 extends ASTBinOp { public ASTPow2() { super(); } @Override protected String opStr(){ return "**" ;} @Override protected ASTOp make() {return new ASTPow2();}
  @Override double op(double d0, double d1) { return Math.pow(d0,d1);}
  @Override String op(String s0, double d1) {throw new IllegalArgumentException("Cannot exponentiate Strings.");}
  @Override String op(double d0, String s1) {throw new IllegalArgumentException("Cannot exponentiate Strings.");}
  @Override String op(String s0, String s1) {throw new IllegalArgumentException("Cannot exponentiate Strings.");}
}
class ASTMod extends ASTBinOp { public ASTMod() { super(); } @Override protected String opStr(){ return "mod"  ;} @Override protected ASTOp make() {return new ASTMod ();}
  @Override double op(double d0, double d1) { return d0%d1;}
  @Override String op(String s0, double d1) {throw new IllegalArgumentException("Cannot mod (%) Strings.");}
  @Override String op(double d0, String s1) {throw new IllegalArgumentException("Cannot exponentiate Strings.");}
  @Override String op(String s0, String s1) {throw new IllegalArgumentException("Cannot exponentiate Strings.");}
}
class ASTLT extends ASTBinOp { public ASTLT() { super(); } @Override protected String opStr(){ return "<"  ;} @Override protected ASTOp make() {return new ASTLT  ();}
  @Override double op(double d0, double d1) { return d0<d1 && !MathUtils.equalsWithinOneSmallUlp(d0,d1)?1:0;}
  @Override String op(String s0, double d1) {throw new IllegalArgumentException("Cannot apply '<' to Strings.");}
  @Override String op(double d0, String s1) {throw new IllegalArgumentException("Cannot apply '<' to Strings.");}
  @Override String op(String s0, String s1) {throw new IllegalArgumentException("Cannot apply '<' to Strings.");}
}
class ASTLE extends ASTBinOp { public ASTLE() { super(); } @Override protected String opStr(){ return "<=" ;} @Override protected ASTOp make() {return new ASTLE  ();}
  @Override double op(double d0, double d1) { return d0<d1 ||  MathUtils.equalsWithinOneSmallUlp(d0,d1)?1:0;}
  @Override String op(String s0, double d1) {throw new IllegalArgumentException("Cannot apply '<=' to Strings.");}
  @Override String op(double d0, String s1) {throw new IllegalArgumentException("Cannot apply '<=' to Strings.");}
  @Override String op(String s0, String s1) {throw new IllegalArgumentException("Cannot apply '<=' to Strings.");}
}
class ASTGT extends ASTBinOp { public ASTGT() { super(); } @Override protected String opStr(){ return ">"  ;} @Override protected ASTOp make() {return new ASTGT  ();}
  @Override double op(double d0, double d1) { return d0>d1 && !MathUtils.equalsWithinOneSmallUlp(d0,d1)?1:0;}
  @Override String op(String s0, double d1) {throw new IllegalArgumentException("Cannot apply '>' to Strings.");}
  @Override String op(double d0, String s1) {throw new IllegalArgumentException("Cannot apply '>' to Strings.");}
  @Override String op(String s0, String s1) {throw new IllegalArgumentException("Cannot apply '>' to Strings.");}
}
class ASTGE extends ASTBinOp { public ASTGE() { super(); } @Override protected String opStr(){ return ">=" ;} @Override protected ASTOp make() {return new ASTGE  ();}
  @Override double op(double d0, double d1) { return d0>d1 ||  MathUtils.equalsWithinOneSmallUlp(d0,d1)?1:0;}
  @Override String op(String s0, double d1) {throw new IllegalArgumentException("Cannot apply '>=' to Strings.");}
  @Override String op(double d0, String s1) {throw new IllegalArgumentException("Cannot apply '>=' to Strings.");}
  @Override String op(String s0, String s1) {throw new IllegalArgumentException("Cannot apply '>=' to Strings.");}
}
class ASTEQ extends ASTBinOp { public ASTEQ() { super(); } @Override protected String opStr(){ return "==" ;} @Override protected ASTOp make() {return new ASTEQ  ();}
  @Override double op(double d0, double d1) { return MathUtils.equalsWithinOneSmallUlp(d0,d1)?1:0;}
  @Override String op(String s0, double d1) { return s0.equals(Double.toString(d1)) ? "1.0" : "0.0"; }
  @Override String op(double d0, String s1) { return (Double.toString(d0)).equals(s1) ? "1.0" : "0.0";}
  @Override String op(String s0, String s1) { return s0.equals(s1) ? "1.0" : "0.0"; }
}
class ASTNE extends ASTBinOp { public ASTNE() { super(); } @Override protected String opStr(){ return "!=" ;} @Override protected ASTOp make() {return new ASTNE  ();}
  @Override double op(double d0, double d1) { return MathUtils.equalsWithinOneSmallUlp(d0,d1)?0:1;}
  @Override String op(String s0, double d1) { return !s0.equals(Double.toString(d1)) ? "1.0" : "0.0"; }
  @Override String op(double d0, String s1) { return !(Double.toString(d0)).equals(s1) ? "1.0" : "0.0";}
  @Override String op(String s0, String s1) { return !s0.equals(s1) ? "1.0" : "0.0"; }
}
class ASTLA extends ASTBinOp { public ASTLA() { super(); } @Override protected String opStr(){ return "&"  ;} @Override protected ASTOp make() {return new ASTLA  ();}
  @Override double op(double d0, double d1) { return (d0!=0 && d1!=0) ? (Double.isNaN(d0) || Double.isNaN(d1)?Double.NaN:1) :0;}
  @Override String op(String s0, double d1) {throw new IllegalArgumentException("Cannot '&' Strings.");}
  @Override String op(double d0, String s1) {throw new IllegalArgumentException("Cannot '&' Strings.");}
  @Override String op(String s0, String s1) {throw new IllegalArgumentException("Cannot '&' Strings.");}
}
class ASTLO extends ASTBinOp { public ASTLO() { super(); } @Override protected String opStr(){ return "|"  ;} @Override protected ASTOp make() {return new ASTLO  ();}
  @Override double op(double d0, double d1) {
    if (d0 == 0 && Double.isNaN(d1)) { return Double.NaN; }
    if (d1 == 0 && Double.isNaN(d0)) { return Double.NaN; }
    if (Double.isNaN(d0) && Double.isNaN(d1)) { return Double.NaN; }
    if (d0 == 0 && d1 == 0) { return 0; }
    return 1;
  }
  @Override String op(String s0, double d1) {throw new IllegalArgumentException("Cannot '|' Strings.");}
  @Override String op(double d0, String s1) {throw new IllegalArgumentException("Cannot '|' Strings.");}
  @Override String op(String s0, String s1) {throw new IllegalArgumentException("Cannot '|' Strings.");}
}

class O extends ASTOp {
  // example (O "a" "sum" "b" "null" (+ %a %b))
  // this is (O _accName _acc _elemName _elem AST)
  String _accName; // the name of the accumulator variable
  String _acc;     // the accumulator
  String _elemName; // the name of the element variable
  String _elem; // the new element, if null, then use the chunk.at(i) value.
  O() { super(null); }
  @Override protected String opStr() { return "O"; }
  @Override protected ASTOp make() { return new O(); }
  protected O parse_impl(Exec E) {
    _accName = E.parseString(E.peekPlus()); E.skipWS();
    _acc     = E.parseString(E.peekPlus()); E.skipWS();
    _elemName= E.parseString(E.peekPlus()); E.skipWS();
    AST elem = E.parse();
    if( elem instanceof ASTNum) _elem=""+((ASTNum)elem)._d;
    else if(elem instanceof ASTString) _elem=((ASTString)elem)._s;
    if( _elem.equals("null")) _elem=null;
    AST ast = E.parse();
    O res = (O)clone();
    res._asts = new AST[]{ast};
    return res;
  }
  @Override protected void apply(Env e) { }
  void exec(NonBlockingHashMap<String,Val> m, Chunk c, int row) {
    Env e = (new Env(null)).capture();
    Val v = m.get(_acc);
    // if v is not null use the type, otherwise use the type of the elem!
    if( v!=null ) e._local.put(_accName, v.type(), v.value());
    int t=Env.NULL;
    if( _elem==null ) {
      if (c.vec().isNumeric())
        e._local.put(_elemName, t=Env.NUM, ""+c.atd(row));
      else if (c.vec().isString())
        e._local.put(_elemName, t=Env.STR, c.atStr(new ValueString(), row).toString());
    } else {
      int type;
      try {
        Double.valueOf(_elem);
        type=Env.NUM;
      } catch(Exception ex) { type=Env.STR; }
      e._local.put(_elemName, t=type, _elem);
    }
    if( v==null )
      e._local.put(_accName, t, t==Env.STR?"":"0");  // currently expects only Strings or Nums...
    _asts[0].treeWalk(e);
    m.put(_acc,e.pop());
  }
  void reduce(NonBlockingHashMap<String,Val> thiz, NonBlockingHashMap<String,Val> that) {
    Env e =(new Env(null)).capture();
    Val l = thiz.get(_acc);
    Val r = that.get(_acc);
    e._local.put(_accName, l.type(),l.value());
    e._local.put(_elemName,r.type(),r.value());
    _asts[0].treeWalk(e);
    thiz.put(_acc, e.pop());
  }
}

class ROp extends ASTOp {
  HashMap<String, O> _ops;
  // parse_impl: (R #N accum1 O accum2 O ...)
  ROp() {super(null); _ops=new HashMap<>(); }
  @Override protected String opStr() { return "R"; }
  @Override protected ASTOp make() { return new ROp(); }
  protected ROp parse_impl(Exec E) {
    double n = ((ASTNum)(E.parse()))._d;
    for(int i=0;i<n;++i) {
      E.skipWS();
      String acc = E.parseString(E.peekPlus()); E.skipWS();
      O o = (O)E.parse();
      _ops.put(acc,o);
    }
    return (ROp)clone();
  }
  void map(NonBlockingHashMap<String,Val> m, Chunk c, int row) {
    for( String s:_ops.keySet() )
      _ops.get(s).exec(m,c,row);
  }
  void reduce(NonBlockingHashMap<String,Val> thiz, NonBlockingHashMap<String,Val> that) {
    for( String s:_ops.keySet())
      _ops.get(s).reduce(thiz,that);
  }
  @Override public AutoBuffer write_impl(AutoBuffer ab) {
    if( _ops==null ) return ab.put4(0);
    ab.put4(_ops.size());
    for( String s:_ops.keySet()) { ab.putStr(s); ab.put(_ops.get(s)); }
    return ab;
  }
  @Override public ROp read_impl(AutoBuffer ab) {
    int len = ab.get4();
    if( len==0 ) return this;
    _ops = new HashMap<>();
    for( int i=0;i<len;++i)
      _ops.put(ab.getStr(), ab.get(O.class));
    return this;
  }
  @Override void exec(Env e) {}
  @Override String value() { return null; }
  @Override int type() { return 0; }
  @Override protected void apply(Env e) {}
}

class COp extends ASTOp {
  COp() {super(null);}
  @Override protected String opStr() { return "C"; }
  @Override protected ASTOp make() { return new COp(); }
  // parse_impl: (C (AST))
  protected COp parse_impl(Exec E) {
    AST ast = E.parse();
    COp res = (COp)clone();
    res._asts = new AST[]{ast};
    return res;
  }
  Val combine(NonBlockingHashMap<String,Val> m) {
    Env e = (new Env(null)).capture();
    for( String s:m.keySet() ) {
      e._local.put(s,m.get(s).type(),m.get(s).value());
    }
    _asts[0].treeWalk(e);
    return e.pop();
  }
  @Override void exec(Env e) {}
  @Override String value() { return null; }
  @Override int type() { return 0; }
  @Override protected void apply(Env e) {}
}

// R like binary operator &&
class ASTAND extends ASTBinOp {
  @Override protected String opStr() { return "&&"; }
  ASTAND( ) { super();}
  @Override double op(double d0, double d1) { throw H2O.fail(); }
  @Override String op(String s0, double d1) {throw new IllegalArgumentException("Cannot '&&' Strings.");}
  @Override String op(double d0, String s1) {throw new IllegalArgumentException("Cannot '&&' Strings.");}
  @Override String op(String s0, String s1) {throw new IllegalArgumentException("Cannot '&&' Strings.");}

  @Override protected ASTOp make() { return new ASTAND(); }
  @Override protected void apply(Env env) {
    double op1 = (env.isNum()) ? env.peekDbl()
            : (env.isAry() ? env.peekAry().vecs()[0].at(0) : Double.NaN);
    env.pop();
    double op2 = (env.isNum()) ? env.peekDbl()
            : (env.isAry() ? env.peekAry().vecs()[0].at(0) : Double.NaN);
    env.pop();

    // Both NAN ? push NaN
    if (Double.isNaN(op1) && Double.isNaN(op2)) {
      env.push(new ValNum(Double.NaN));
      return;
    }

    // Either 0 ? push False
    if (op1 == 0 || op2 == 0) {
      env.push(new ValNum(0.0));
      return;
    }

    // Either NA ? push NA (no need to worry about 0s, taken care of in case above)
    if (Double.isNaN(op1) || Double.isNaN(op2)) {
      env.push(new ValNum(Double.NaN));
      return;
    }

    // Otherwise, push True
    env.push(new ValNum(1.0));
  }
}

// R like binary operator ||
class ASTOR extends ASTBinOp {
  @Override protected String opStr() { return "||"; }
  ASTOR( ) { super(); }
  @Override double op(double d0, double d1) { throw H2O.fail(); }
  @Override String op(String s0, double d1) {throw new IllegalArgumentException("Cannot '||' Strings.");}
  @Override String op(double d0, String s1) {throw new IllegalArgumentException("Cannot '||' Strings.");}
  @Override String op(String s0, String s1) {throw new IllegalArgumentException("Cannot '||' Strings.");}

  @Override protected ASTOp make() { return new ASTOR(); }
  @Override protected void apply(Env env) {
    double op1 = (env.isNum()) ? env.peekDbl()
            : (env.isAry() ? env.peekAry().vecs()[0].at(0) : Double.NaN);
    // op1 is NaN ? push NaN
    if (Double.isNaN(op1)) {
      env.poppush(2, new ValNum(Double.NaN));
      return;
    }
    double op2 = !Double.isNaN(op1) && op1!=0 ? 1 : (env.isNum()) ? env.peekDbl()
                    : (env.isAry()) ? env.peekAry().vecs()[0].at(0) : Double.NaN;

    // op2 is NaN ? push NaN
    if (Double.isNaN(op2)) {
      env.poppush(2, new ValNum(op2));
      return;
    }

    // both 0 ? push False
    if (op1 == 0 && op2 == 0) {
      env.poppush(2, new ValNum(0.0));
      return;
    }

    // else push True
    env.poppush(2, new ValNum(1.0));
  }
}

/**
* R 'ls' command.
*
* This method is purely for the console right now.  Print stuff into the string buffer.
* JSON response is not configured at all.
*/
class ASTLs extends ASTOp {
  ASTLs() { super(new String[]{"ls"}); }
  @Override protected String opStr() { return "ls"; }
  @Override protected ASTOp make() {return new ASTLs();}
  protected ASTLs parse_impl(Exec E) { return (ASTLs) clone(); }
  @Override protected void apply(Env env) {
    ArrayList<String> domain = new ArrayList<>();
    Futures fs = new Futures();
    AppendableVec av = new AppendableVec(Vec.VectorGroup.VG_LEN1.addVec());
    NewChunk keys = new NewChunk(av,0);
    int r = 0;
    for( Key key : KeySnapshot.globalSnapshot().keys()) {
      keys.addEnum(r++);
      domain.add(key.toString());
    }
    keys.close(fs);
    Vec c0 = av.close(fs);   // c0 is the row index vec
    fs.blockForPending();
    String[] key_domain = new String[domain.size()];
    for (int i = 0; i < key_domain.length; ++i) key_domain[i] = domain.get(i);
    c0.setDomain(key_domain);
    env.pushAry(new Frame(Key.make("h2o_ls"), new String[]{"key"}, new Vec[]{c0}));
  }

  private double getSize(Key k) {
    return (double)(((Frame) k.get()).byteSize());
//    if (k.isChunkKey()) return (double)((Chunk)DKV.get(k).get()).byteSize();
//    if (k.isVec()) return (double)((Vec)DKV.get(k).get()).rollupStats()._size;
//    return Double.NaN;
  }
}


// WIP

class ASTXorSum extends ASTReducerOp {
  ASTXorSum() {super(0); }
  @Override protected String opStr(){ return "xorsum";}
  @Override protected ASTOp make() {return new ASTXorSum();}
  @Override protected double op(double d0, double d1) {
    long d0Bits = Double.doubleToLongBits(d0);
    long d1Bits = Double.doubleToLongBits(d1);
    long xorsumBits = d0Bits ^ d1Bits;
    // just need to not get inf or nan. If we zero the upper 4 bits, we won't
    final long ZERO_SOME_SIGN_EXP = 0x0fffffffffffffffL;
    xorsumBits = xorsumBits & ZERO_SOME_SIGN_EXP;
    return Double.longBitsToDouble(xorsumBits);
  }
  @Override double[] map(Env env, double[] in, double[] out, AST[] args) {
    if (out == null || out.length < 1) out = new double[1];
    long xorsumBits = 0;
    long vBits;
    // for dp ieee 754 , sign and exp are the high 12 bits
    // We don't want infinity or nan, because h2o will return a string.
    double xorsum = 0;
    for (double v : in) {
      vBits = Double.doubleToLongBits(v);
      xorsumBits = xorsumBits ^ vBits;
    }
    // just need to not get inf or nan. If we zero the upper 4 bits, we won't
    final long ZERO_SOME_SIGN_EXP = 0x0fffffffffffffffL;
    xorsumBits = xorsumBits & ZERO_SOME_SIGN_EXP;
    xorsum = Double.longBitsToDouble(xorsumBits);
    out[0] = xorsum;
    return out;
  }
}

class ASTMMult extends ASTOp {
  ASTMMult() { super(VARS2);}

  protected ASTMMult parse_impl(Exec E) {
    AST l = E.parse();
    if (l instanceof ASTId) l = Env.staticLookup((ASTId)l);
    AST r = E.parse();
    if (r instanceof ASTId) r = Env.staticLookup((ASTId)r);
    ASTMMult res = new ASTMMult();
    res._asts = new AST[]{l,r};
    return res;
  }
  @Override
  protected String opStr() { return "x";}

  @Override
  protected ASTOp make() { return new ASTMMult();}

  @Override
  protected void apply(Env env) {
    env.poppush(2, new ValFrame(DMatrix.mmul(env.peekAryAt(-0), env.peekAryAt(-1))));
  }
}

class ASTTranspose extends ASTOp {
  ASTTranspose() { super(VARS1);}

  protected ASTTranspose parse_impl(Exec E) {
    AST arg = E.parse();
    if (arg instanceof ASTId) arg = Env.staticLookup((ASTId)arg);
    ASTTranspose res = new ASTTranspose();
    res._asts = new AST[]{arg};
    return res;
  }
  @Override
  protected String opStr() { return "t";}

  @Override
  protected ASTOp make() { return new ASTTranspose();}

  @Override
  protected void apply(Env env) {
    env.push(new ValFrame(DMatrix.transpose(env.popAry())));
  }
}


// Legacy Items: On the chopping block


// Brute force implementation of matrix multiply
//class ASTMMult extends ASTOp {
//  @Override protected String opStr() { return "%*%"; }
//  ASTMMult( ) {
//    super(new String[]{"", "x", "y"},
//            new Type[]{Type.ARY,Type.ARY,Type.ARY},
//            OPF_PREFIX,
//            OPP_MUL,
//            OPA_RIGHT);
//  }
//  @Override protected ASTOp make() { return new ASTMMult(); }
//  @Override protected void apply(Env env, int argcnt, ASTApply apply) {
//    env.poppush(3,new Matrix(env.ary(-2)).mult(env.ary(-1)),null);
//  }
//}
//
//// Brute force implementation of matrix transpose
//class ASTMTrans extends ASTOp {
//  @Override protected String opStr() { return "t"; }
//  ASTMTrans( ) {
//    super(new String[]{"", "x"},
//            new Type[]{Type.ARY,Type.dblary()},
//            OPF_PREFIX,
//            OPP_PREFIX,
//            OPA_RIGHT);
//  }
//  @Override protected ASTOp make() { return new ASTMTrans(); }
//  @Override protected void apply(Env env, int argcnt, ASTApply apply) {
//    if(!env.isAry(-1)) {
//      Key k = new Vec.VectorGroup().addVec();
//      Futures fs = new Futures();
//      AppendableVec avec = new AppendableVec(k);
//      NewChunk chunk = new NewChunk(avec, 0);
//      chunk.addNum(env.dbl(-1));
//      chunk.close(0, fs);
//      Vec vec = avec.close(fs);
//      fs.blockForPending();
//      vec.domain = null;
//      Frame fr = new Frame(new String[] {"C1"}, new Vec[] {vec});
//      env.poppush(2,new Matrix(fr).trans(),null);
//    } else
//      env.poppush(2,new Matrix(env.ary(-1)).trans(),null);
//  }
//}


//class ASTPrint extends ASTOp {
//  static Type[] newsig() {
//    Type t1 = Type.unbound();
//    return new Type[]{t1, t1, Type.varargs(Type.unbound())};
//  }
//  ASTPrint() { super(new String[]{"print", "x", "y..."},
//          newsig(),
//          OPF_PREFIX,
//          OPP_PREFIX,OPA_RIGHT); }
//  @Override protected String opStr() { return "print"; }
//  @Override protected ASTOp make() {return new ASTPrint();}
//  @Override protected void apply(Env env, int argcnt, ASTApply apply) {
//    for( int i=1; i<argcnt; i++ ) {
//      if( env.isAry(i-argcnt) ) {
//        env._sb.append(env.ary(i-argcnt).toStringAll());
//      } else {
//        env._sb.append(env.toString(env._sp+i-argcnt,true));
//      }
//    }
//    env.pop(argcnt-2);          // Pop most args
//    env.pop_into_stk(-2);       // Pop off fcn, returning 1st arg
//  }
//}



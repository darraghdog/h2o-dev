package water.rapids;

import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.*;
import java.util.concurrent.atomic.AtomicInteger;
import jsr166y.CountedCompleter;
import org.apache.commons.math3.special.Gamma;
import org.apache.commons.math3.util.FastMath;
import org.joda.time.MutableDateTime;
import sun.misc.Unsafe;
import water.*;
import water.exceptions.H2OIllegalArgumentException;
import water.fvec.*;
import water.nbhm.NonBlockingHashMap;
import water.nbhm.NonBlockingHashSet;
import water.nbhm.UtilUnsafe;
import water.parser.ValueString;
import water.util.ArrayUtils;
import water.util.Log;

public abstract class ASTUniPrefixOp extends ASTUniOp {
  public ASTUniPrefixOp( ) { super(); }
  public ASTUniPrefixOp( String[] vars) { super(vars); }
}

class ASTNot  extends ASTUniPrefixOp { public ASTNot()  { super(); } @Override protected String opStr(){ return "!";} @Override protected ASTOp make() {return new ASTNot(); } @Override double op(double d) { if (Double.isNaN(d)) return Double.NaN; return d==0?1:0; } }
class ASTCos  extends ASTUniPrefixOp { @Override protected String opStr(){ return "cos";  } @Override protected ASTOp make() {return new ASTCos ();} @Override double op(double d) { return Math.cos(d);}}
class ASTSin  extends ASTUniPrefixOp { @Override protected String opStr(){ return "sin";  } @Override protected ASTOp make() {return new ASTSin ();} @Override double op(double d) { return Math.sin(d);}}
class ASTTan  extends ASTUniPrefixOp { @Override protected String opStr(){ return "tan";  } @Override protected ASTOp make() {return new ASTTan ();} @Override double op(double d) { return Math.tan(d);}}
class ASTACos extends ASTUniPrefixOp { @Override protected String opStr(){ return "acos"; } @Override protected ASTOp make() {return new ASTACos();} @Override double op(double d) { return Math.acos(d);}}
class ASTASin extends ASTUniPrefixOp { @Override protected String opStr(){ return "asin"; } @Override protected ASTOp make() {return new ASTASin();} @Override double op(double d) { return Math.asin(d);}}
class ASTATan extends ASTUniPrefixOp { @Override protected String opStr(){ return "atan"; } @Override protected ASTOp make() {return new ASTATan();} @Override double op(double d) { return Math.atan(d);}}
class ASTCosh extends ASTUniPrefixOp { @Override protected String opStr(){ return "cosh"; } @Override protected ASTOp make() {return new ASTCosh ();} @Override double op(double d) { return Math.cosh(d);}}
class ASTSinh extends ASTUniPrefixOp { @Override protected String opStr(){ return "sinh"; } @Override protected ASTOp make() {return new ASTSinh ();} @Override double op(double d) { return Math.sinh(d);}}
class ASTTanh extends ASTUniPrefixOp { @Override protected String opStr(){ return "tanh"; } @Override protected ASTOp make() {return new ASTTanh ();} @Override double op(double d) { return Math.tanh(d);}}
class ASTACosh extends ASTUniPrefixOp { @Override protected String opStr(){ return "acosh"; } @Override protected ASTOp make() {return new ASTACosh ();} @Override double op(double d) { return FastMath.acosh(d);}}
class ASTASinh extends ASTUniPrefixOp { @Override protected String opStr(){ return "asinh"; } @Override protected ASTOp make() {return new ASTASinh ();} @Override double op(double d) { return FastMath.asinh(d);}}
class ASTATanh extends ASTUniPrefixOp { @Override protected String opStr(){ return "atanh"; } @Override protected ASTOp make() {return new ASTATanh ();} @Override double op(double d) { return FastMath.atanh(d);}}
class ASTCosPi extends ASTUniPrefixOp { @Override protected String opStr(){ return "cospi"; } @Override protected ASTOp make() {return new ASTCosPi ();} @Override double op(double d) { return Math.cos(Math.PI*d);}}
class ASTSinPi extends ASTUniPrefixOp { @Override protected String opStr(){ return "sinpi"; } @Override protected ASTOp make() {return new ASTSinPi ();} @Override double op(double d) { return Math.sin(Math.PI*d);}}
class ASTTanPi extends ASTUniPrefixOp { @Override protected String opStr(){ return "tanpi"; } @Override protected ASTOp make() {return new ASTTanPi ();} @Override double op(double d) { return Math.tan(Math.PI*d);}}
class ASTAbs  extends ASTUniPrefixOp { @Override protected String opStr(){ return "abs";  } @Override protected ASTOp make() {return new ASTAbs ();} @Override double op(double d) { return Math.abs(d);}}
class ASTSgn  extends ASTUniPrefixOp { @Override protected String opStr(){ return "sign" ; } @Override protected ASTOp make() {return new ASTSgn ();} @Override double op(double d) { return Math.signum(d);}}
class ASTSqrt extends ASTUniPrefixOp { @Override protected String opStr(){ return "sqrt"; } @Override protected ASTOp make() {return new ASTSqrt();} @Override double op(double d) { return Math.sqrt(d);}}
class ASTTrun extends ASTUniPrefixOp { @Override protected String opStr(){ return "trunc"; } @Override protected ASTOp make() {return new ASTTrun();} @Override double op(double d) { return d>=0?Math.floor(d):Math.ceil(d);}}
class ASTCeil extends ASTUniPrefixOp { @Override protected String opStr(){ return "ceiling"; } @Override protected ASTOp make() {return new ASTCeil();} @Override double op(double d) { return Math.ceil(d);}}
class ASTFlr  extends ASTUniPrefixOp { @Override protected String opStr(){ return "floor";} @Override protected ASTOp make() {return new ASTFlr ();} @Override double op(double d) { return Math.floor(d);}}
class ASTLog  extends ASTUniPrefixOp { @Override protected String opStr(){ return "log";  } @Override protected ASTOp make() {return new ASTLog ();} @Override double op(double d) { return Math.log(d);}}
class ASTLog10  extends ASTUniPrefixOp { @Override protected String opStr(){ return "log10";  } @Override protected ASTOp make() {return new ASTLog10 ();} @Override double op(double d) { return Math.log10(d);}}
class ASTLog2  extends ASTUniPrefixOp { @Override protected String opStr(){ return "log2";  } @Override protected ASTOp make() {return new ASTLog2 ();} @Override double op(double d) { return Math.log(d)/Math.log(2);}}
class ASTLog1p  extends ASTUniPrefixOp { @Override protected String opStr(){ return "log1p";  } @Override protected ASTOp make() {return new ASTLog1p ();} @Override double op(double d) { return Math.log1p(d);}}
class ASTExp  extends ASTUniPrefixOp { @Override protected String opStr(){ return "exp";  } @Override protected ASTOp make() {return new ASTExp ();} @Override double op(double d) { return Math.exp(d);}}
class ASTExpm1  extends ASTUniPrefixOp { @Override protected String opStr(){ return "expm1";  } @Override protected ASTOp make() {return new ASTExpm1 ();} @Override double op(double d) { return Math.expm1(d);}}
class ASTGamma  extends ASTUniPrefixOp { @Override protected String opStr(){ return "gamma";  } @Override protected ASTOp make() {return new ASTGamma ();} @Override double op(double d) {  return Gamma.gamma(d);}}
class ASTLGamma extends ASTUniPrefixOp { @Override protected String opStr(){ return "lgamma"; } @Override protected ASTOp make() {return new ASTLGamma ();} @Override double op(double d) { return Gamma.logGamma(d);}}
class ASTDiGamma  extends ASTUniPrefixOp { @Override protected String opStr(){ return "digamma";  } @Override protected ASTOp make() {return new ASTDiGamma ();} @Override double op(double d) {  return Gamma.digamma(d);}}
class ASTTriGamma  extends ASTUniPrefixOp { @Override protected String opStr(){ return "trigamma";  } @Override protected ASTOp make() {return new ASTTriGamma ();} @Override double op(double d) {  return Gamma.trigamma(d);}}

class ASTIsNA extends ASTUniPrefixOp { @Override protected String opStr(){ return "is.na";} @Override protected ASTOp make() { return new ASTIsNA();} @Override double op(double d) { return Double.isNaN(d)?1:0;}
  @Override protected void apply(Env env) {
    // Expect we can broadcast across all functions as needed.
    if( env.isNum() ) { env.push(new ValNum(op(env.popDbl()))); return; }
    Frame fr = env.popAry();
    Frame fr2 = new MRTask() {
      @Override public void map( Chunk chks[], NewChunk nchks[] ) {
        for( int i=0; i<nchks.length; i++ ) {
          NewChunk n = nchks[i];
          Chunk c = chks[i];
          int rlen = c._len;
          for( int r=0; r<rlen; r++ )
            n.addNum( c.isNA(r) ? 1 : 0);
        }
      }
    }.doAll(fr.numCols(),fr).outputFrame(Key.make(), fr._names, null);
    env.pushAry(fr2);
  }
}

class ASTRound extends ASTUniPrefixOp {
  int _digits = 0;
  @Override protected String opStr() { return "round"; }
  ASTRound() { super(new String[]{"round", "x", "digits"}); }
  @Override
  protected ASTRound parse_impl(Exec E) {
    // Get the ary
    if (!E.hasNext()) throw new IllegalArgumentException("End of input unexpected. Badly formed AST.");
    AST ary = E.parse();
    if (ary instanceof ASTId) ary = Env.staticLookup((ASTId)ary);
    // Get the digits
    if (!(E.skipWS().hasNext())) throw new IllegalArgumentException("End of input unexpected. Badly formed AST.");
    try {
      _digits = (int) ((ASTNum) (E.parse())).dbl();
    } catch (ClassCastException e) {
      e.printStackTrace();
      throw new IllegalArgumentException("Expected a number for `digits` argument.");
    }
    ASTRound res = (ASTRound) clone();
    res._asts = new AST[]{ary};
    return res;
  }
  @Override protected ASTOp make() { return new ASTRound(); }
  @Override protected void apply(Env env) {
    final int digits = _digits;
    if(env.isAry()) {
      Frame fr = env.popAry();
      for(int i = 0; i < fr.vecs().length; i++) {
        if(fr.vecs()[i].isEnum())
          throw new IllegalArgumentException("Non-numeric column " + String.valueOf(i+1) + " in data frame");
      }
      Frame fr2 = new MRTask() {
        @Override public void map(Chunk chks[], NewChunk nchks[]) {
          for(int i = 0; i < nchks.length; i++) {
            NewChunk n = nchks[i];
            Chunk c = chks[i];
            int rlen = c._len;
            for(int r = 0; r < rlen; r++)
              n.addNum(roundDigits(c.atd(r),digits));
          }
        }
      }.doAll(fr.numCols(),fr).outputFrame(fr.names(),fr.domains());
      env.pushAry(fr2);
    }
    else
      env.push(new ValNum(roundDigits(env.popDbl(), digits)));
  }

  // e.g.: floor(2.676*100 + 0.5) / 100 => 2.68
  static double roundDigits(double x, int digits) {
    if(Double.isNaN(x)) return x;
    double sgn = x < 0 ? -1 : 1;
    x = Math.abs(x);
    double power_of_10 = (int)Math.pow(10, digits);
    return sgn*(digits == 0
            // go to the even digit
            ? (x % 1 >= 0.5 && !(Math.floor(x)%2==0))
              ? Math.ceil(x)
              : Math.floor(x)
            : Math.floor(x * power_of_10 + 0.5) / power_of_10);
  }
}

class ASTSignif extends ASTUniPrefixOp {
  int _digits = 6;  // R default
  @Override protected String opStr() { return "signif"; }
  ASTSignif() { super(new String[]{"signif", "x", "digits"}); }
  @Override
  protected ASTSignif parse_impl(Exec E) {
    // Get the ary
    AST ary = E.parse();
    if (ary instanceof ASTId) ary = Env.staticLookup((ASTId)ary);
    // Get the digits
    try {
      _digits = (int) ((ASTNum) (E.parse())).dbl();
    } catch (ClassCastException e) {
      e.printStackTrace();
      throw new IllegalArgumentException("Expected a double for `digits` argument.");
    }
    ASTSignif res = (ASTSignif) clone();
    res._asts = new AST[]{ary};
    return res;
  }
  @Override protected ASTOp make() { return new ASTSignif(); }
  @Override protected void apply(Env env) {
    final int digits = _digits;
    if(digits < 0)
      throw new IllegalArgumentException("Error in signif: argument digits must be a non-negative integer");

    if(env.isAry()) {
      Frame fr = env.popAry();
      for(int i = 0; i < fr.vecs().length; i++) {
        if(fr.vecs()[i].isEnum())
          throw new IllegalArgumentException("Non-numeric column " + String.valueOf(i+1) + " in data frame");
      }
      Frame fr2 = new MRTask() {
        @Override public void map(Chunk chks[], NewChunk nchks[]) {
          for(int i = 0; i < nchks.length; i++) {
            NewChunk n = nchks[i];
            Chunk c = chks[i];
            int rlen = c._len;
            for(int r = 0; r < rlen; r++)
              n.addNum(signifDigits(c.atd(r),digits));
          }
        }
      }.doAll(fr.numCols(),fr).outputFrame(fr.names(),fr.domains());
      env.pushAry(fr2);
    }
    else
      env.push(new ValNum(signifDigits(env.popDbl(), digits)));
  }
  static double signifDigits(double x, int digits) {
    if(Double.isNaN(x)) return x;
    BigDecimal bd = new BigDecimal(x);
    bd = bd.round(new MathContext(digits, RoundingMode.HALF_EVEN));
    return bd.doubleValue();
  }
}

class ASTNrow extends ASTUniPrefixOp {
  public ASTNrow() { super(VARS1); }
  @Override protected String opStr() { return "nrow"; }
  @Override protected ASTOp make() {return new ASTNrow();}
  @Override protected void apply(Env env) {
    Frame fr = env.popAry();
    double d = fr.numRows();
    env.push(new ValNum(d));
  }
}

class ASTNcol extends ASTUniPrefixOp {
  ASTNcol() { super(VARS1); }
  @Override protected String opStr() { return "ncol"; }
  @Override protected ASTOp make() {return new ASTNcol();}
  @Override protected void apply(Env env) {
    Frame fr = env.popAry();
    double d = fr.numCols();
    env.push(new ValNum(d));
  }
}

class ASTLength extends ASTUniPrefixOp {
  ASTLength() { super(VARS1); }
  @Override protected String opStr() { return "length"; }
  @Override protected ASTOp make() { return new ASTLength(); }
  @Override protected void apply(Env env) {
    Frame fr = env.popAry();
    double d = fr.numCols() == 1 ? fr.numRows() : fr.numCols();
//    env.cleanup(fr);
//    env.clean();
    env.push(new ValNum(d));
  }
}

class ASTIsFactor extends ASTUniPrefixOp {
  ASTIsFactor() { super(VARS1); }
  @Override protected String opStr() { return "is.factor"; }
  @Override protected ASTOp make() {return new ASTIsFactor();}
  @Override protected void apply(Env env) {
    Frame fr = env.popAry();
    String res = "FALSE";
    if (fr.numCols() != 1) throw new IllegalArgumentException("is.factor applies to a single column.");
    if (fr.anyVec().isEnum()) res = "TRUE";
    env.push(new ValStr(res));
  }
}

// Added to facilitate Runit testing
class ASTAnyFactor extends ASTUniPrefixOp {
  ASTAnyFactor() { super(VARS1);}
  @Override protected String opStr() { return "any.factor"; }
  @Override protected ASTOp make() {return new ASTAnyFactor();}
  @Override protected void apply(Env env) {
    Frame fr = env.popAry();
    String res = "FALSE";
    for (int i = 0; i < fr.vecs().length; ++i)
      if (fr.vecs()[i].isEnum()) { res = "TRUE"; break; }
    env.push(new ValStr(res));
  }
}

class ASTCanBeCoercedToLogical extends ASTUniPrefixOp {
  ASTCanBeCoercedToLogical() { super(VARS1); }
  @Override protected String opStr() { return "canBeCoercedToLogical"; }
  @Override protected ASTOp make() {return new ASTCanBeCoercedToLogical();}
  @Override protected void apply(Env env) {
    Frame fr = env.popAry();
    String res = "FALSE";
    Vec[] v = fr.vecs();
    for (Vec aV : v)
      if (aV.isInt())
        if (aV.min() == 0 && aV.max() == 1) { res = "TRUE"; break; }
    env.push(new ValStr(res));
  }
}

class ASTAnyNA extends ASTUniPrefixOp {
  ASTAnyNA() { super(VARS1); }
  @Override protected String opStr() { return "any.na"; }
  @Override protected ASTOp make() {return new ASTAnyNA();}
  @Override protected void apply(Env env) {
    Frame fr = env.popAry();
    String res = "FALSE";
    for (int i = 0; i < fr.vecs().length; ++i)
      if (fr.vecs()[i].naCnt() > 0) { res = "TRUE"; break; }
    env.push(new ValStr(res));
  }
}

//class ASTIsTRUE extends ASTUniPrefixOp {
//  ASTIsTRUE() {super(VARS1,new Type[]{Type.DBL,Type.unbound()});}
//  @Override protected String opStr() { return "isTRUE"; }
//  @Override protected ASTOp make() {return new ASTIsTRUE();}  // to make sure fcn get bound at each new context
//  @Override protected void apply(Env env, int argcnt, ASTApply apply) {
//    double res = env.isDbl() && env.popDbl()==1.0 ? 1:0;
//    env.pop();
//    env.poppush(res);
//  }
//}

class ASTScale extends ASTUniPrefixOp {
  boolean _center;
  double[] _centers;
  boolean _scale;
  double[] _scales;
  ASTScale() { super(new String[]{"ary", "center", "scale"});}
  @Override protected String opStr() { return "scale"; }
  @Override protected ASTOp make() {return new ASTScale();}
  @Override
  protected ASTScale parse_impl(Exec E) {
    AST ary = E.parse();
    if (ary instanceof ASTId) ary = Env.staticLookup((ASTId)ary);
    parseArg(E, true);  // centers parse
    parseArg(E, false); // scales parse
    ASTScale res = (ASTScale) clone();
    res._asts = new AST[]{ary};
    return res;
  }
  private void parseArg(Exec E, boolean center) {
    if (center) {
      if (!E.skipWS().hasNext()) throw new IllegalArgumentException("End of input unexpected. Badly formed AST.");
      String[] centers = E.peek() == '{' ? E.xpeek('{').parseString('}').split(";") : null;
      if (centers == null) {
        // means `center` is boolean
        AST a;
        try {
          a = E._env.lookup((ASTId) E.skipWS().parse());
        } catch (ClassCastException e) {
          e.printStackTrace();
          throw new IllegalArgumentException("Expected to get an ASTId. Badly formed AST.");
        }
        try {
          _center = ((ASTNum) a).dbl() == 1;
          _centers = null;
        } catch (ClassCastException e) {
          e.printStackTrace();
          throw new IllegalArgumentException("Expected to get a number for the `center` argument.");
        }
      } else {
        for (int i = 0; i < centers.length; ++i) centers[i] = centers[i].replace("\"", "").replace("\'", "");
        _centers = new double[centers.length];
        for (int i = 0; i < centers.length; ++i) _centers[i] = Double.valueOf(centers[i]);
      }
    } else {
      if (!E.skipWS().hasNext()) throw new IllegalArgumentException("End of input unexpected. Badly formed AST.");
      String[] centers = E.peek() == '{' ? E.xpeek('{').parseString('}').split(";") : null;
      if (centers == null) {
        // means `scale` is boolean
        AST a;
        try {
          a = E._env.lookup((ASTId) E.skipWS().parse());
        } catch (ClassCastException e) {
          e.printStackTrace();
          throw new IllegalArgumentException("Expected to get an ASTId. Badly formed AST.");
        }
        try {
          _scale = ((ASTNum) a).dbl() == 1;
          _scales = null;
        } catch (ClassCastException e) {
          e.printStackTrace();
          throw new IllegalArgumentException("Expected to get a number for the `scale` argument.");
        }
      } else {
        for (int i = 0; i < centers.length; ++i) centers[i] = centers[i].replace("\"", "").replace("\'", "");
        _scales = new double[centers.length];
        for (int i = 0; i < centers.length; ++i) _scales[i] = Double.valueOf(centers[i]);
      }
    }
  }

  @Override protected void apply(Env env) {
    Frame fr = env.popAry();
    for (int i = 0; i < fr.numCols(); ++i) if (fr.vecs()[i].isEnum()) throw new IllegalArgumentException(("All columns must be numeric."));
    if (!(_centers == null) && _centers.length != fr.numCols()) throw new IllegalArgumentException("`centers` must be logical or have length equal to the number of columns in the dataset.");
    if (!(_scales  == null) && _scales.length  != fr.numCols()) throw new IllegalArgumentException("`scales` must be logical or have length equal to the number of columns in the dataset.");
    final boolean use_mean = _centers == null && _center;
    final double[] centers = _centers;
    final boolean use_sig  = _scales == null && _scale;
    final boolean use_rms  = !use_mean && _scale;
    final double[] scales  = _scales;
    if (!_center && !_scale && (_centers == null) && (_scales == null)) {
      //nothing to do, return the frame as is
      env.pushAry(fr);
      return;
    }

    boolean doCenter = use_mean || _centers != null;
    boolean doScale  = use_sig || use_rms || _scales != null;

    Frame centered = new Frame(fr);
    if (doCenter) {
      centered = new MRTask() {
        @Override
        public void map(Chunk[] cs, NewChunk[] ncs) {
          int rows = cs[0]._len;
          int cols = cs.length;
          for (int r = 0; r < rows; ++r)
            for (int c = 0; c < cols; ++c) {
              double numer = cs[c].atd(r) - (use_mean
                      ? cs[c].vec().mean()
                      : centers == null ? 0 : centers[c]);
              ncs[c].addNum(numer);
            }
        }
      }.doAll(fr.numCols(), fr).outputFrame(fr.names(), fr.domains());
    }

    double[] rms_vals = null;
    if (use_rms) {
      rms_vals = new double[fr.numCols()];
      double nrows = fr.numRows();
      for (int i = 0; i < rms_vals.length; ++i) {
        Vec v = centered.vecs()[i];
        ASTVar.CovarTask t = new ASTVar.CovarTask(0,0).doAll(new Frame(v,v));
        rms_vals[i] = Math.sqrt(t._ss / (nrows - 1));
      }
    }
    final double[] rms = rms_vals;

    Frame scaled = new Frame(centered);
    if (doScale) {
      scaled = new MRTask() {
        @Override
        public void map(Chunk[] cs, NewChunk[] ncs) {
          int rows = cs[0]._len;
          int cols = cs.length;
          for (int r = 0; r < rows; ++r)
            for (int c = 0; c < cols; ++c) {
              double denom = cs[c].atd(r) / (use_rms
                      ? rms[c] : use_sig ? cs[c].vec().sigma()
                      : scales == null ? 1 : scales[c]);
              ncs[c].addNum(denom);
            }
        }
      }.doAll(centered.numCols(), centered).outputFrame(centered.names(), centered.domains());
    }
    env.pushAry(scaled);
  }
}

abstract class ASTTimeOp extends ASTUniPrefixOp {
  ASTTimeOp() { super(VARS1); }
  // Override for e.g. month and day-of-week
  protected String[][] factors() { return null; }
  @Override
  protected ASTTimeOp parse_impl(Exec E) {
    AST arg = E.parse();
    if (arg instanceof ASTId) arg = Env.staticLookup((ASTId)arg);
    ASTTimeOp res = (ASTTimeOp) clone();
    res._asts = new AST[]{arg};
    return res;
  }

  abstract long op( MutableDateTime dt );

  @Override protected void apply(Env env) {
    // Single instance of MDT for the single call
    if( !env.isAry() ) {        // Single point
      double d = env.peekDbl();
      if( !Double.isNaN(d) ) d = op(new MutableDateTime((long)d));
      env.poppush(1, new ValNum(d));
      return;
    }
    // Whole column call
    Frame fr = env.peekAry();
    final ASTTimeOp uni = this;
    Frame fr2 = new MRTask() {
      @Override public void map( Chunk chks[], NewChunk nchks[] ) {
        MutableDateTime dt = new MutableDateTime(0);
        for( int i=0; i<nchks.length; i++ ) {
          NewChunk n =nchks[i];
          Chunk c = chks[i];
          int rlen = c._len;
          for( int r=0; r<rlen; r++ ) {
            double d = c.atd(r);
            if( !Double.isNaN(d) ) {
              dt.setMillis((long)d);
              d = uni.op(dt);
            }
            n.addNum(d);
          }
        }
      }
    }.doAll(fr.numCols(),fr).outputFrame(fr._names, factors());
    env.poppush(1, new ValFrame(fr2));
  }
}
//
class ASTYear  extends ASTTimeOp { @Override protected String opStr(){ return "year" ; } @Override protected ASTOp make() {return new ASTYear  ();} @Override long op(MutableDateTime dt) { return dt.getYear();}}
class ASTDay   extends ASTTimeOp { @Override protected String opStr(){ return "day"  ; } @Override protected ASTOp make() {return new ASTDay   ();} @Override long op(MutableDateTime dt) { return dt.getDayOfMonth();}}
class ASTHour  extends ASTTimeOp { @Override protected String opStr(){ return "hour" ; } @Override protected ASTOp make() {return new ASTHour  ();} @Override long op(MutableDateTime dt) { return dt.getHourOfDay();}}
class ASTMinute extends ASTTimeOp { @Override protected String opStr(){return "minute";} @Override protected ASTOp make() {return new ASTMinute();} @Override long op(MutableDateTime dt) { return dt.getMinuteOfHour();}}
class ASTSecond extends ASTTimeOp { @Override protected String opStr(){return "second";} @Override protected ASTOp make() {return new ASTSecond();} @Override long op(MutableDateTime dt) { return dt.getSecondOfMinute();}}
class ASTMillis extends ASTTimeOp { @Override protected String opStr(){return "millis";} @Override protected ASTOp make() {return new ASTMillis();} @Override long op(MutableDateTime dt) { return dt.getMillisOfSecond();}}
class ASTMonth extends ASTTimeOp { 
  static private final String[][] FACTORS = new String[][]{{"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"}};
  @Override protected String[][] factors() { return FACTORS; }
  @Override protected String opStr(){ return "month"; } 
  @Override protected ASTOp make() {return new ASTMonth ();} 
  @Override long op(MutableDateTime dt) { return dt.getMonthOfYear()-1;}
}
class ASTDayOfWeek extends ASTTimeOp { 
  static private final String[][] FACTORS = new String[][]{{"Mon","Tue","Wed","Thu","Fri","Sat","Sun"}}; // Order comes from Joda
  @Override protected String[][] factors() { return FACTORS; }
  @Override protected String opStr(){ return "dayOfWeek"; } 
  @Override protected ASTOp make() {return new ASTDayOfWeek();} 
  @Override long op(MutableDateTime dt) { return dt.getDayOfWeek()-1;}
}

// Convert year, month, day, hour, minute, sec, msec to Unix epoch time
class ASTMktime extends ASTUniPrefixOp {
  ASTMktime() { super(new String[]{"","year","month","day","hour","minute","second","msec"}); }
  @Override protected String opStr() { return "mktime"; }
  @Override protected ASTMktime make() {return new ASTMktime();} 
  @Override
  protected ASTMktime parse_impl(Exec E) {
    AST yr = E.parse();  if( yr instanceof ASTId) yr = Env.staticLookup((ASTId)yr);
    AST mo = E.parse();  if( mo instanceof ASTId) mo = Env.staticLookup((ASTId)mo);
    AST dy = E.parse();  if( dy instanceof ASTId) dy = Env.staticLookup((ASTId)dy);
    AST hr = E.parse();  if( hr instanceof ASTId) hr = Env.staticLookup((ASTId)hr);
    AST mi = E.parse();  if( mi instanceof ASTId) mi = Env.staticLookup((ASTId)mi);
    AST se = E.parse();  if( se instanceof ASTId) se = Env.staticLookup((ASTId)se);
    AST ms = E.parse();  if( ms instanceof ASTId) ms = Env.staticLookup((ASTId)ms);
    ASTMktime res = (ASTMktime) clone();
    res._asts = new AST[]{yr,mo,dy,hr,mi,se,ms};
    return res;
  }

  @Override protected void apply(Env env) {
    // Seven args, all required.  See if any are arrays.
    Frame fs[] = new Frame[7];
    int   is[] = new int  [7];
    Frame x = null;             // Sample frame (for auto-expanding constants)
    for( int i=0; i<7; i++ )
      if( env.peekType()==Env.ARY ) fs[i] = x = env.popAry();
      else                          is[i] =(int)env.popDbl();

    if( x==null ) {                            // Single point
      long msec = new MutableDateTime(is[6],   // year   
                                      is[5]+1, // month  
                                      is[4]+1, // day    
                                      is[3],   // hour   
                                      is[2],   // minute 
                                      is[1],   // second 
                                      is[0])   // msec   
        .getMillis();
      env.poppush(1, new ValNum(msec));
      return;
    }

    // Make constant Vecs for the constant args.  Commonly, they'll all be zero
    Vec vecs[] = new Vec[7];
    for( int i=0; i<7; i++ ) {
      if( fs[i] == null ) {
        vecs[i] = x.anyVec().makeCon(is[i]);
      } else {
        if( fs[i].numCols() != 1 ) throw new IllegalArgumentException("Expect single column");
        vecs[i] = fs[i].anyVec();
      }
    }

    // Convert whole column to epoch msec
    Frame fr2 = new MRTask() {
      @Override public void map( Chunk chks[], NewChunk nchks[] ) {
        MutableDateTime dt = new MutableDateTime(0);
        NewChunk n = nchks[0];
        int rlen = chks[0]._len;
        for( int r=0; r<rlen; r++ ) {
          dt.setDateTime((int)chks[6].at8(r),  // year   
                         (int)chks[5].at8(r)+1,// month  
                         (int)chks[4].at8(r)+1,// day    
                         (int)chks[3].at8(r),  // hour   
                         (int)chks[2].at8(r),  // minute 
                         (int)chks[1].at8(r),  // second 
                         (int)chks[0].at8(r)); // msec   
          n.addNum(dt.getMillis());
        }
      }
      }.doAll(1,vecs).outputFrame(new String[]{"msec"},null);
    env.poppush(1, new ValFrame(fr2));
  }
}

// operate on a single vec
// reduce the Vec
class ASTFoldCombine extends ASTUniPrefixOp {
  // (RC (R ...) (C ...) vec)
  private ROp _red;     // operates on a single value
  private COp _combine; // what to do with the _accum map
  ASTFoldCombine() { super(null); }
  @Override protected String opStr() { return "RC"; }
  @Override protected ASTOp make() { return new ASTFoldCombine(); }
  protected ASTFoldCombine parse_impl(Exec E) {
    _red = (ROp)E.parse();
    _combine = (COp)E.parse();
    AST ary =  E.parse();
    ASTFoldCombine res = (ASTFoldCombine) clone();
    res._asts = new AST[]{ary};
    return res;
  }
  @Override protected void apply(Env e) {
    Frame f = e.popAry();
    if( f.numCols() != 1 )
      throw new IllegalArgumentException("Expected one column, got "+f.numCols());

    // apply _red across the frame f and fetch out _accum
    NonBlockingHashMap<String,Val> accum = new RTask(_red).doAll(f.anyVec())._accum;

    // then apply the _combine operator on the accum...
    e.push(_combine.combine(accum));
  }

  private static class RTask extends MRTask<RTask> {
    NonBlockingHashMap<String,Val> _accum;
    private ROp _red;
    RTask(ROp red) { _red=red; }
    @Override public void setupLocal() {_accum=new NonBlockingHashMap<>();}
    @Override public void map(Chunk cs) {
      for( int i=0;i<cs._len;++i)
        _red.map(_accum,cs,i);
    }
    @Override public void reduce(RTask t) { _red.reduce(_accum,t._accum); }
    @Override public AutoBuffer write_impl( AutoBuffer ab ) {
      ab.put(_red);
      if( _accum == null ) return ab.put4(0);
      ab.put4(_accum.size());
      for( String s:_accum.keySet() ) {
        ab.putStr(s);
        ab.put(_accum.get(s));
      }
      return ab;
    }
    @Override public RTask read_impl(AutoBuffer ab) {
      _red = ab.get(ROp.class);
      int len = ab.get4();
      if( len == 0 ) return this;
      _accum = new NonBlockingHashMap<>();
      for( int i=0;i<len;++i )
        _accum.put(ab.getStr(), ab.get(Val.class));
      return this;
    }
  }
}

class ASTRbind extends ASTUniPrefixOp {
  protected static int argcnt;
  @Override protected String opStr() { return "rbind"; }
  public ASTRbind() { super(new String[]{"rbind", "ary","..."}); }
  @Override protected ASTOp make() { return new ASTRbind(); }
  protected ASTRbind parse_impl(Exec E) {
    ArrayList<AST> dblarys = new ArrayList<>();
    AST ary = E.parse();
    if (ary instanceof ASTId) ary = Env.staticLookup((ASTId)ary);
    dblarys.add(ary);
    AST a;
    boolean broke = false;
    while (E.skipWS().hasNext()) {
      a = E.parse();
      if (a instanceof ASTId) {
        AST ast = E._env.lookup((ASTId)a);
        if (ast instanceof ASTFrame) { dblarys.add(a); }
        else {broke = true; break; } // if not a frame then break here since we are done parsing Frame args
      }
      else if (a instanceof ASTFrame || a instanceof ASTSlice || a instanceof ASTBinOp || a instanceof ASTUniOp || a instanceof ASTReducerOp) { // basically anything that returns a Frame...
        dblarys.add(a);
      }
      else { broke = true; break; }
    }
    if (broke) E.rewind();
    Collections.reverse(dblarys);
    ASTRbind res = (ASTRbind) clone();
    res._asts = dblarys.toArray(new AST[argcnt=dblarys.size()]);
    return res;
  }

  private String get_type(byte t) {
    switch(t) {
      case Vec.T_ENUM: return "factor";
      case Vec.T_NUM:  return "numeric";
      case Vec.T_STR:  return "String";
      case Vec.T_TIME: return "time";
      case Vec.T_UUID: return "UUID";
      default: return "bad";
    }
  }

  private static class RbindMRTask extends MRTask<RbindMRTask> {
    private final int[] _emap;
    private final int _chunkOffset;
    private final Vec _v;
    RbindMRTask(H2O.H2OCountedCompleter hc, int[] emap, Vec v, int offset) { super(hc); _emap = emap; _v = v; _chunkOffset = offset;}

    @Override public void map(Chunk cs) {
      int idx = _chunkOffset+cs.cidx();
      Key ckey = Vec.chunkKey(_v._key, idx);
      if (_emap != null) {
        assert !cs.hasFloat(): "Input chunk ("+cs.getClass()+") has float, but is expected to be enum";
        NewChunk nc = new NewChunk(_v, idx);
        // loop over rows and update ints for new domain mapping according to vecs[c].domain()
        for (int r=0;r < cs._len;++r) {
          if (cs.isNA(r)) nc.addNA();
          else nc.addNum(_emap[(int)cs.at8(r)], 0);
        }
        nc.close(_fs);
      } else {
        Chunk oc = (Chunk)cs.clone();
        oc.setStart(-1);
        oc.setVec(null);
        oc.setBytes(cs.getBytes().clone()); // needless replication of the data, can do ref counting on byte[] _mem
        DKV.put(ckey, oc, _fs, true);
      }
    }
  }

  private static class RbindTask extends H2O.H2OCountedCompleter<RbindTask> {
    final transient Vec[] _vecs;
    final Vec _v;
    final long[] _espc;
    String[] _dom;

    RbindTask(H2O.H2OCountedCompleter cc, Vec[] vecs, Vec v, long[] espc) { super(cc); _vecs = vecs; _v = v; _espc = espc; }

    private static Map<Integer, String> invert(Map<String, Integer> map) {
      Map<Integer, String> inv = new HashMap<>();
      for (Map.Entry<String, Integer> e : map.entrySet()) {
        inv.put(e.getValue(), e.getKey());
      }
      return inv;
    }

    @Override protected void compute2() {
      addToPendingCount(_vecs.length-1);
      boolean isEnum = _vecs[0].domain() != null;
      int[][] emaps  = new int[_vecs.length][];

      if (isEnum) {
        // loop to create BIG domain
        HashMap<String, Integer> dmap = new HashMap<>(); // probably should allocate something that's big enough (i.e. 2*biggest_domain)
        int c = 0;
        for (int i = 0; i < _vecs.length; ++i) {
          emaps[i] = new int[_vecs[i].domain().length];
          for (int j = 0; j < emaps[i].length; ++j)
            if (!dmap.containsKey(_vecs[i].domain()[j]))
              dmap.put(_vecs[i].domain()[j], emaps[i][j]=c++);
            else emaps[i][j] = dmap.get(_vecs[i].domain()[j]);
        }
        _dom = new String[dmap.size()];
        HashMap<Integer, String> inv = (HashMap<Integer, String>) invert(dmap);
        for (int s = 0; s < _dom.length; ++s) _dom[s] = inv.get(s);
      }
      int offset=0;
      for (int i=0; i<_vecs.length; ++i) {
        new RbindMRTask(this, emaps[i], _v, offset).asyncExec(_vecs[i]);
        offset += _vecs[i].nChunks();
      }
    }

    @Override public void onCompletion(CountedCompleter cc) {
        _v.setDomain(_dom);
        DKV.put(_v);
    }
  }

  private static class ParallelRbinds extends H2O.H2OCountedCompleter{

    private final Env _env;
    private final int _argcnt;
    private final AtomicInteger _ctr;
    private int _maxP = 100;

    private long[] _espc;
    private Vec[] _vecs;
    ParallelRbinds(Env e, int argcnt) { _env = e; _argcnt = argcnt; _ctr = new AtomicInteger(_maxP-1); }  //TODO pass maxP to constructor

    @Override protected void compute2() {
      addToPendingCount(_env.peekAry().numCols()-1);
      int nchks=0;
      for (int i =0; i < _argcnt; ++i)
        nchks+=_env.peekAryAt(-i).anyVec().nChunks();

      _espc = new long[nchks+1];
      int coffset = _env.peekAry().anyVec().nChunks();
      long[] first_espc = _env.peekAry().anyVec().get_espc();
      System.arraycopy(first_espc, 0, _espc, 0, first_espc.length);
      for (int i=1; i< _argcnt; ++i) {
        long roffset = _espc[coffset];
        long[] espc = _env.peekAryAt(-i).anyVec().get_espc();
        int j = 1;
        for (; j < espc.length; j++)
          _espc[coffset + j] = roffset+ espc[j];
        coffset += _env.peekAryAt(-i).anyVec().nChunks();
      }

      Key[] keys = _env.peekAry().anyVec().group().addVecs(_env.peekAry().numCols());
      _vecs = new Vec[keys.length];
      for (int i=0; i<_vecs.length; ++i)
        _vecs[i] = new Vec( keys[i], _espc, null, _env.peekAry().vec(i).get_type());

      for (int i=0; i < Math.min(_maxP, _vecs.length); ++i) forkVecTask(i);
    }

    private void forkVecTask(final int i) {
      Vec[] vecs = new Vec[_argcnt];
      for (int j= 0; j < _argcnt; ++j)
        vecs[j] = _env.peekAryAt(-j).vec(i);
      new RbindTask(new Callback(), vecs, _vecs[i], _espc).fork();
    }

    private class Callback extends H2O.H2OCallback {
      public Callback(){super(ParallelRbinds.this);}
      @Override public void callback(H2O.H2OCountedCompleter h2OCountedCompleter) {
        int i = _ctr.incrementAndGet();
        if(i < _vecs.length)
          forkVecTask(i);
      }
    }
  }

  @Override protected void apply(Env env) {
    // quick check to make sure rbind is feasible
    if (argcnt == 1) { return; } // leave stack as is

    Frame f1 = env.peekAry();
    // do error checking and compute new offsets in tandem
    for (int i = 1; i < argcnt; ++i) {
      Frame t = env.peekAryAt(-i);

      // check columns match
      if (t.numCols() != f1.numCols())
        throw new IllegalArgumentException("Column mismatch! Expected " + f1.numCols() + " but frame has " + t.numCols());

      // check column types
      for (int c = 0; c < f1.numCols(); ++c) {
        if (f1.vec(c).get_type() != t.vec(c).get_type())
          throw new IllegalArgumentException("Column type mismatch! Expected type " + get_type(f1.vec(c).get_type()) + " but vec has type " + get_type(t.vec(c).get_type()));
      }
    }
    ParallelRbinds t;
    H2O.submitTask(t =new ParallelRbinds(env, argcnt)).join();
    env.poppush(argcnt, new ValFrame(new Frame(f1.names(), t._vecs)));
  }
}

class ASTCbind extends ASTUniPrefixOp {
  protected static int argcnt;
  @Override protected String opStr() { return "cbind"; }
  public ASTCbind() { super(new String[]{"cbind","ary", "..."}); }
  @Override protected ASTOp make() {return new ASTCbind();}
  protected ASTCbind parse_impl(Exec E) {
    ArrayList<AST> dblarys = new ArrayList<>();
    AST ary = E.parse();
    if (ary instanceof ASTId) ary = Env.staticLookup((ASTId)ary);
    dblarys.add(ary);
    AST a;
    boolean broke = false;
    while (E.skipWS().hasNext()) {
      a = E.parse();
      if (a instanceof ASTId) {
        AST ast = E._env.lookup((ASTId)a);
        if (ast instanceof ASTFrame) { dblarys.add(a); }
        else {broke = true; break; } // if not a frame then break here since we are done parsing Frame args
      }
      else if (a instanceof ASTFrame || a instanceof ASTSlice || a instanceof ASTBinOp || a instanceof ASTUniPrefixOp || a instanceof ASTUniOp || a instanceof ASTReducerOp) { // basically anything that returns a Frame...
        dblarys.add(a);
      }
      else { broke = true; break; }
    }
    if (broke) E.rewind();
    ASTCbind res = (ASTCbind) clone();
    AST[] arys = new AST[argcnt=dblarys.size()];
    for (int i = 0; i < dblarys.size(); i++) arys[i] = dblarys.get(i);
    res._asts = arys;
    return res;
  }
  @Override protected void apply(Env env) {
    // Validate the input frames
    Vec vmax = null;
    for(int i = 0; i < argcnt; i++) {
      Frame t = env.peekAryAt(-i);
      if(vmax == null) vmax = t.vecs()[0];
      else if(t.numRows() != vmax.length())
        // R pads shorter cols to match max rows by cycling/repeating, but we won't support that
        throw new IllegalArgumentException("Row mismatch! Expected " + String.valueOf(vmax.length()) + " but frame has " + String.valueOf(t.numRows()));
    }

    // loop over frames and combine
    Frame fr = new Frame(new String[0],new Vec[0]);
    for(int i = 0; i < argcnt; i++) {
      Frame f = env.peekAryAt(i-argcnt+1);  // Reverse order off stack
      Frame ff = f.deepSlice(null,null);  // deep copy the frame, R semantics...
      Frame new_frame = fr.makeCompatible(ff);
      if (f.numCols() == 1) fr.add(f.names()[0], new_frame.anyVec());
      else fr.add(new_frame);
    }
    env.pop(argcnt);

    env.pushAry(fr);
  }
}

class ASTRename extends ASTUniPrefixOp {
  protected static String _newname;
  @Override protected String opStr() { return "rename"; }
  ASTRename() { super(new String[] {"", "ary", "new_name"}); }
  @Override protected ASTOp make() { return new ASTRename(); }
  protected ASTRename parse_impl(Exec E) {
    AST ary = E.parse();
    if (ary instanceof ASTId) ary = Env.staticLookup((ASTId)ary);
    _newname = ((ASTString)E.parse())._s;
    ASTRename res = (ASTRename) clone();
    res._asts = new AST[]{ary};
    return res;
  }

  @Override protected void apply(Env e) {
    Frame fr = e.popAry();
    Frame fr2 = fr.deepCopy(_newname);
    DKV.put(fr2._key, fr2);
    e.pushAry(fr2);
  }
}

class ASTSetLevel extends ASTUniPrefixOp {
  private String _lvl;
  ASTSetLevel() { super(new String[]{"setLevel", "x", "level"});}
  @Override protected String opStr() { return "setLevel"; }
  @Override protected ASTOp make() { return new ASTSetLevel(); }
  protected ASTSetLevel parse_impl(Exec E) {
    AST ary = E.parse();
    if( ary instanceof ASTId ) ary = Env.staticLookup((ASTId)ary);
    _lvl = ((ASTString)E.parse())._s;
    ASTSetLevel res = (ASTSetLevel) clone();
    res._asts = new AST[]{ary};
    return res;
  }
  @Override protected void apply(Env env) {
    Frame fr = env.peekAry();
    if (fr.numCols() != 1) throw new IllegalArgumentException("`setLevel` works on a single column at a time.");
    String[] doms = fr.anyVec().domain().clone();
    if( doms == null )
      throw new IllegalArgumentException("Cannot set the level on a non-factor column!");
    final int idx = Arrays.asList(doms).indexOf(_lvl);
    if (idx == -1)
      throw new IllegalArgumentException("Did not find level `" + _lvl + "` in the column.");

    Frame fr2 = new MRTask() {
      @Override public void map(Chunk c, NewChunk nc) {
        for (int i=0;i<c._len;++i)
          nc.addNum(idx);
      }
    }.doAll(1, fr.anyVec()).outputFrame(null, fr.names(), fr.domains());
    env.poppush(1, new ValFrame(fr2));
  }
}

class ASTMatch extends ASTUniPrefixOp {
  protected static double _nomatch;
  protected static String[] _matches;
  @Override protected String opStr() { return "match"; }
  ASTMatch() { super( new String[]{"", "ary", "table", "nomatch", "incomparables"}); }
  @Override protected ASTOp make() { return new ASTMatch(); }
  protected ASTMatch parse_impl(Exec E) {
    // First parse out the `ary` arg
    AST ary = E.parse();
    if (ary instanceof ASTId) ary = Env.staticLookup((ASTId)ary);
    // The `table` arg
    if (!E.skipWS().hasNext()) throw new IllegalArgumentException("End of input unexpected. Badly formed AST.");
    _matches = E.peek() == '{' ? E.xpeek('{').parseString('}').split(";") : new String[]{E.parseString(E.peekPlus())};
    // cleanup _matches
    for (int i = 0; i < _matches.length; ++i) _matches[i] = _matches[i].replace("\"", "").replace("\'", "");
    // `nomatch` is just a number in case no match
    try {
      ASTNum nomatch = (ASTNum) E.skipWS().parse();
      _nomatch = nomatch.dbl();
    } catch(ClassCastException e) {
      e.printStackTrace();
      throw new IllegalArgumentException("Argument `nomatch` expected a number.");
    }
    // drop the incomparables arg for now ...
    AST incomp = E.skipWS().parse();
    ASTMatch res = (ASTMatch) clone();
    res._asts = new AST[]{ary};
    return res;
  }

  @Override protected void apply(Env e) {
    Frame fr = e.popAry();
    if (fr.numCols() != 1) throw new IllegalArgumentException("can only match on a single categorical column.");
    if (!fr.anyVec().isEnum()) throw new IllegalArgumentException("can only match on a single categorical column.");
    Key tmp = Key.make();
    final String[] matches = _matches;
    Frame rez = new MRTask() {
      private int in(String s) { return Arrays.asList(matches).contains(s) ? 1 : 0; }
      @Override public void map(Chunk c, NewChunk n) {
        int rows = c._len;
        for (int r = 0; r < rows; ++r) n.addNum(in(c.vec().domain()[(int)c.at8(r)]));
      }
    }.doAll(1, fr.anyVec()).outputFrame(tmp, null, null);
    e.pushAry(rez);
  }

}

// Similar to R's seq_len
class ASTSeqLen extends ASTUniPrefixOp {
  protected static double _length;
  @Override protected String opStr() { return "seq_len"; }
  ASTSeqLen( ) { super(new String[]{"seq_len", "n"}); }
  @Override protected ASTOp make() { return new ASTSeqLen(); }
  @Override
  protected ASTSeqLen parse_impl(Exec E) {
    try {
      _length = E.nextDbl();
    } catch(ClassCastException e) {
      e.printStackTrace();
      throw new IllegalArgumentException("Argument `n` expected to be a number.");
    }
    ASTSeqLen res = (ASTSeqLen) clone();
    res._asts = new AST[]{};
    return res;
  }

  @Override protected void apply(Env env) {
    int len = (int) Math.ceil(_length);
    if (len <= 0)
      throw new IllegalArgumentException("Error in seq_len(" +len+"): argument must be coercible to positive integer");
    Frame fr = new Frame(new String[]{"c"}, new Vec[]{Vec.makeSeq(len)});
    env.pushAry(fr);
  }
}

// Same logic as R's generic seq method
class ASTSeq extends ASTUniPrefixOp {
  protected static double _from;
  protected static double _to;
  protected static double _by;

  @Override protected String opStr() { return "seq"; }
  ASTSeq() { super(new String[]{"seq", "from", "to", "by"}); }
  @Override protected ASTOp make() { return new ASTSeq(); }
  @Override
  protected ASTSeq parse_impl(Exec E) {
    // *NOTE*: This function creates a frame, there is no input frame!

    // Get the from
    try {
      if (!E.skipWS().hasNext()) throw new IllegalArgumentException("End of input unexpected. Badly formed AST. Missing `from` argument.");
      _from = E.nextDbl();
    } catch (ClassCastException e) {
      e.printStackTrace();
      throw new IllegalArgumentException("Argument `from` expected to be a number.");
    }
    // Get the to
    try {
      if (!E.skipWS().hasNext()) throw new IllegalArgumentException("End of input unexpected. Badly formed AST. Missing `to` argument.");
      _to = E.nextDbl();
    } catch (ClassCastException e) {
      e.printStackTrace();
      throw new IllegalArgumentException("Argument `to` expected to be a number.");
    }
    // Get the by
    try {
      if (!E.skipWS().hasNext()) throw new IllegalArgumentException("End of input unexpected. Badly formed AST. Missing `by` argument.");
      _by = E.nextDbl();
    } catch (ClassCastException e) {
      e.printStackTrace();
      throw new IllegalArgumentException("Argument `by` expected to be a number.");
    }
    // Finish the rest
    ASTSeq res = (ASTSeq) clone();
    res._asts = new AST[]{}; // in reverse order so they appear correctly on the stack.
    return res;
  }

  @Override protected void apply(Env env) {
    double delta = _to - _from;
    if(delta == 0 && _to == 0)
      env.push(new ValNum(_to));
    else {
      double n = delta/_by;
      if(n < 0)
        throw new IllegalArgumentException("wrong sign in 'by' argument");
      else if(n > Double.MAX_VALUE)
        throw new IllegalArgumentException("'by' argument is much too small");

      double dd = Math.abs(delta)/Math.max(Math.abs(_from), Math.abs(_to));
      if(dd < 100*Double.MIN_VALUE)
        env.push(new ValNum(_from));
      else {
        Futures fs = new Futures();
        AppendableVec av = new AppendableVec(Vec.newKey());
        NewChunk nc = new NewChunk(av, 0);
        int len = (int)n + 1;
        for (int r = 0; r < len; r++) nc.addNum(_from + r*_by);
        // May need to adjust values = by > 0 ? min(values, to) : max(values, to)
        nc.close(0, fs);
        Vec vec = av.close(fs);
        fs.blockForPending();
        Frame fr = new Frame(new String[]{"C1"}, new Vec[]{vec});
        env.pushAry(fr);
      }
    }
  }
}

class ASTRepLen extends ASTUniPrefixOp {
  protected static double _length;
  @Override protected String opStr() { return "rep_len"; }
  public ASTRepLen() { super(new String[]{"rep_len", "x", "length.out"}); }
  @Override protected ASTOp make() { return new ASTRepLen(); }
  protected ASTRepLen parse_impl(Exec E) {
    AST ary = E.parse();
    if (ary instanceof ASTId) { ary = Env.staticLookup((ASTId)ary); }
    try {
      _length = E.nextDbl();
    } catch(ClassCastException e) {
      e.printStackTrace();
      throw new IllegalArgumentException("Argument `length` expected to be a number.");
    }
    ASTRepLen res = (ASTRepLen) clone();
    res._asts = new AST[]{ary};
    return res;
  }
  @Override protected void apply(Env env) {

    // two cases if x is a frame: x is a single vec, x is a list of vecs
    if (env.isAry()) {
      final Frame fr = env.popAry();
      if (fr.numCols() == 1) {

        // In this case, create a new vec of length _length using the elements of x
        Vec v = Vec.makeRepSeq((long)_length, fr.numRows());  // vec of "indices" corresponding to rows in x
        new MRTask() {
          @Override public void map(Chunk c) {
            for (int i = 0; i < c._len; ++i)
              c.set(i, fr.anyVec().at((long) c.atd(i)));
          }
        }.doAll(v);
        v.setDomain(fr.anyVec().domain());
        Frame f = new Frame(new String[]{"C1"}, new Vec[]{v});
        env.pushAry(f);

      } else {

        // In this case, create a new frame with numCols = _length using the vecs of fr
        // this is equivalent to what R does, but by additionally calling "as.data.frame"
        String[] col_names = new String[(int)_length];
        for (int i = 0; i < col_names.length; ++i) col_names[i] = "C" + (i+1);
        Frame f = new Frame(col_names, new Vec[(int)_length]);
        for (int i = 0; i < f.numCols(); ++i)
          f.add(Frame.defaultColName(f.numCols()), fr.vec( i % fr.numCols() ));
        env.pushAry(f);
      }
    }

    // x is a number or a string
    else {
      int len = (int)_length;
      if(len <= 0)
        throw new IllegalArgumentException("Error in rep_len: argument length.out must be coercible to a positive integer");
      if (env.isStr()) {
        // make a constant enum vec with domain[] = []{env.popStr()}
        Frame fr = new Frame(new String[]{"C1"}, new Vec[]{Vec.makeCon(0, len)});
        fr.anyVec().setDomain(new String[]{env.popStr()});
        env.pushAry(fr);
      } else if (env.isNum()) {
        Frame fr = new Frame(new String[]{"C1"}, new Vec[]{Vec.makeCon(env.popDbl(), len)});
        env.pushAry(fr);
      } else throw new IllegalArgumentException("Unkown input. Type: "+env.peekType() + " Stack: " + env.toString());
    }
  }
}

class ASTSetColNames extends ASTUniPrefixOp {
  protected static long[] _idxs;
  protected static String[] _names;
  @Override protected String opStr() { return "colnames="; }
  public ASTSetColNames() { super(new String[]{}); }
  @Override protected ASTSetColNames make() { return new ASTSetColNames(); }

  // AST: (colnames<- $ary {indices} {names})
  // example:  (colnames<- $iris {#3;#5} {new_name1;new_name2})
  // also acceptable: (colnames<- $iris (: #3 #5) {new_name1;new_name2})
  @Override
  protected ASTSetColNames parse_impl(Exec E) {
    // frame we're changing column names of
    AST ary = E.parse();
    if (ary instanceof ASTId) ary = Env.staticLookup((ASTId)ary);
    // col ids: can be a {#;#;#} or (: # #)
    AST a = E.skipWS().parse();
    if (a instanceof ASTSpan || a instanceof ASTSeries) {
      _idxs = (a instanceof ASTSpan) ? ((ASTSpan) a).toArray() : ((ASTSeries) a).toArray();
      Arrays.sort(_idxs);
    } else if (a instanceof ASTNum) {
      _idxs = new long[]{(long)((ASTNum) a).dbl()};
    } else throw new IllegalArgumentException("Bad AST: Expected a span, series, or number for the column indices.");

    // col names must either be an ASTSeries or a single string
    _names = E.skipWS().peek() == '{' ? E.xpeek('{').parseString('}').replaceAll("\"","").split(";") : new String[]{E.parseString(E.peekPlus())};
    if (_names.length != _idxs.length)
      throw new IllegalArgumentException("Mismatch! Number of columns to change ("+(_idxs.length)+") does not match number of names given ("+(_names.length)+").");

    ASTSetColNames res = (ASTSetColNames)clone();
    res._asts = new AST[]{ary};
    return res;
  }

  @Override protected void apply(Env env) {
    Frame f = env.popAry();
    for (int i=0; i < _names.length; ++i)
      f._names[(int)_idxs[i]] = _names[i];
    if (f._key != null && DKV.get(f._key) != null) DKV.put(f);
    env.pushAry(f);
  }
}

class ASTRunif extends ASTUniPrefixOp {
  protected static long   _seed;
  @Override protected String opStr() { return "h2o.runif"; }
  public ASTRunif() { super(new String[]{"h2o.runif","dbls","seed"}); }
  @Override protected ASTOp make() {return new ASTRunif();}
  @Override
  protected ASTRunif parse_impl(Exec E) {
    // peel off the ary
    AST ary = E.parse();
    if (ary instanceof ASTId) ary = Env.staticLookup((ASTId)ary);
    // parse the seed
    try {
      _seed = (long) E.nextDbl();
    } catch (ClassCastException e) {
      e.printStackTrace();
      throw new IllegalArgumentException("Argument `seed` expected to be a number.");
    }
    ASTRunif res = (ASTRunif) clone();
    res._asts = new AST[]{ary};
    return res;
  }

  @Override protected void apply(Env env) {
    final long seed = _seed == -1 ? (new Random().nextLong()) : _seed;
    Vec rnd = env.popAry().anyVec().makeRand(seed);
    Frame f = new Frame(new String[]{"rnd"}, new Vec[]{rnd});
    env.pushAry(f);
  }
}

class ASTSdev extends ASTUniPrefixOp {
  boolean _narm = false;
  public ASTSdev() { super(new String[]{"sd", "ary", "na.rm"}); }
  @Override protected String opStr() { return "sd"; }
  @Override protected ASTOp make() { return new ASTSdev(); }
  @Override
  protected ASTSdev parse_impl(Exec E) {
    // Get the ary
    AST ary = E.parse();
    if (ary instanceof ASTId) ary = Env.staticLookup((ASTId)ary);
    // Get the na.rm
    AST a = E._env.lookup((ASTId)E.skipWS().parse());
    _narm = ((ASTNum)a).dbl() == 1;
    ASTSdev res = (ASTSdev) clone();
    res._asts = new AST[]{ary}; // in reverse order so they appear correctly on the stack.
    return res;
  }
  @Override protected void apply(Env env) {
    if (env.isNum()) {
      env.poppush(1, new ValNum(Double.NaN));
    } else {
      Frame fr = env.peekAry();
      if (fr.vecs().length > 1)
        throw new IllegalArgumentException("sd does not apply to multiple cols.");
      if (fr.vecs()[0].isEnum())
        throw new IllegalArgumentException("sd only applies to numeric vector.");

      double sig = Math.sqrt(ASTVar.getVar(fr.anyVec(), _narm));
      env.poppush(1, new ValNum(sig));
    }
  }
}

class ASTVar extends ASTUniPrefixOp {
  boolean _narm = false;
  boolean _ynull = false;
  public ASTVar() { super(new String[]{"var", "ary", "y", "na.rm", "use"}); } // the order Vals appear on the stack
  @Override protected String opStr() { return "var"; }
  @Override protected ASTOp make() { return new ASTVar(); }
  @Override
  protected ASTVar parse_impl(Exec E) {
    // Get the ary
    AST ary = E.parse();
    if (ary instanceof ASTId) ary = Env.staticLookup((ASTId)ary);
    // Get the trim
    AST y = E.skipWS().parse();
    if (y instanceof ASTString && ((ASTString)y)._s.equals("null")) {_ynull = true; y = ary; }
    else if (y instanceof ASTId) y = Env.staticLookup((ASTId)y);
    // Get the na.rm
    AST a = E._env.lookup((ASTId)E.skipWS().parse());
    try {
      _narm = ((ASTNum) a).dbl() == 1;
    } catch (ClassCastException e) {
      e.printStackTrace();
      throw new IllegalArgumentException("Argument `na.rm` expected to be a number.");
    }
    // Get the `use`
    ASTString use;
    try {
      use = (ASTString) E.skipWS().parse();
    } catch (ClassCastException e) {
      e.printStackTrace();
      throw new IllegalArgumentException("Argument `use` expected to be a string.");
    }
    // Finish the rest
    ASTVar res = (ASTVar) clone();
    res._asts = new AST[]{use,y,ary}; // in reverse order so they appear correctly on the stack.
    return res;
  }

  @Override protected void apply(Env env) {
    if (env.isNum()) {
      env.pop();
      env.push(new ValNum(Double.NaN));
    } else {
      Frame fr = env.peekAry();                   // number of rows
      Frame y = ((ValFrame) env.peekAt(-1))._fr;  // number of columns
      String use = ((ValStr) env.peekAt(-2))._s;  // what to do w/ NAs: "everything","all.obs","complete.obs","na.or.complete","pairwise.complete.obs"
//      String[] rownames = fr.names();  TODO: Propagate rownames?
      String[] colnames = y.names();

      if (fr.numRows() != y.numRows())
        throw new IllegalArgumentException("In var(): incompatible dimensions. Frames must have the same number of rows.");

      if (use.equals("everything")) _narm = false;
      if (use.equals("complete.obs")) _narm = true;
      if (use.equals("all.obs")) _narm = false;

      if (fr.numRows() == 1) {
        double xmean=0; double ymean=0; double divideby = fr.numCols()-1; double ss=0;
        for (Vec v : fr.vecs()) xmean+= v.at(0);
        for (Vec v : y.vecs())  ymean+= v.at(0);
        xmean /= (divideby+1); ymean /= (divideby+1);

        if (Double.isNaN(xmean) || Double.isNaN(ymean)) { ss = Double.NaN; }
        else {
          for (int r = 0; r <= divideby; ++r) {
            ss += (fr.vecs()[r].at(0) - xmean) * (y.vecs()[r].at(0) - ymean);
          }
        }
        env.poppush(3, new ValNum(Double.isNaN(ss) ? ss : ss/divideby));

      } else {

        final double[/*cols*/][/*rows*/] covars = new double[y.numCols()][fr.numCols()];
        final CovarTask tsks[][] = new CovarTask[y.numCols()][fr.numCols()];
        final Frame frs[][] = new Frame[y.numCols()][fr.numCols()];
        final double xmeans[] = new double[fr.numCols()];
        final double ymeans[] = new double[y.numCols()];
        for (int r = 0; r < fr.numCols(); ++r) xmeans[r] = getMean(fr.vecs()[r], _narm, use);
        for (int c = 0; c < y.numCols(); ++c) ymeans[c] = getMean(y.vecs()[c], _narm, use);
        for (int c = 0; c < y.numCols(); ++c) {
          for (int r = 0; r < fr.numCols(); ++r) {
            frs[c][r] = new Frame(y.vecs()[c], fr.vecs()[r]);
            tsks[c][r] = new CovarTask(ymeans[c], xmeans[r]).doAll(frs[c][r]);
          }
        }
        for (int c = 0; c < y.numCols(); c++)
          for (int r = 0; r < fr.numCols(); r++) {
            covars[c][r] = tsks[c][r].getResult()._ss / (fr.numRows() - 1);
            frs[c][r] = null;
          }

        // Just push the scalar if input is a single col
        if (covars.length == 1 && covars[0].length == 1) env.poppush(3, new ValNum(covars[0][0]));
        else {
          // Build output vecs for var-cov matrix
          Key keys[] = Vec.VectorGroup.VG_LEN1.addVecs(covars.length);
          Vec[] vecs = new Vec[covars.length];
          Futures fs = new Futures();
          for (int i = 0; i < covars.length; i++) {
            AppendableVec v = new AppendableVec(keys[i]);
            NewChunk c = new NewChunk(v, 0);
            for (int j = 0; j < covars[0].length; j++) c.addNum(covars[i][j]);
            c.close(0, fs);
            vecs[i] = v.close(fs);
          }
          fs.blockForPending();
          env.poppush(3, new ValFrame(new Frame(colnames, vecs)));
        }
      }
    }
  }

  static double getMean(Vec v, boolean narm, String use) {
    ASTMean.MeanNARMTask t = new ASTMean.MeanNARMTask(narm).doAll(v);
    if (t._rowcnt == 0 || Double.isNaN(t._sum)) {
      if (use.equals("all.obs")) throw new IllegalArgumentException("use = \"all.obs\" with missing observations.");
      return Double.NaN;
    }
    return t._sum / t._rowcnt;
  }

  static double getVar(Vec v, boolean narm) {
    double m = getMean( v, narm, "");
    CovarTask t = new CovarTask(m,m).doAll(new Frame(v, v));
    return t._ss / (v.length() - 1);
  }

  static class CovarTask extends MRTask<CovarTask> {
    double _ss;
    double _xmean;
    double _ymean;
    CovarTask(double xmean, double ymean) { _xmean = xmean; _ymean = ymean; }
    @Override public void map(Chunk[] cs) {
      int len = cs[0]._len;
      Chunk x = cs[0];
      Chunk y = cs[1];
      if (Double.isNaN(_xmean) || Double.isNaN(_ymean)) { _ss = Double.NaN; return; }
      for (int r = 0; r < len; ++r) {
        _ss += (x.atd(r) - _xmean) * (y.atd(r) - _ymean);
      }
    }
    @Override public void reduce(CovarTask tsk) { _ss += tsk._ss; }
  }
}

class ASTMean extends ASTUniPrefixOp {
  double  _trim = 0;
  boolean _narm = false;
  public ASTMean() { super(new String[]{"mean", "ary", "trim", "na.rm"}); }
  @Override protected String opStr() { return "mean"; }
  @Override protected ASTOp make() { return new ASTMean(); }
  @Override
  protected ASTMean parse_impl(Exec E) {
    // Get the ary
    AST ary = E.parse();
    if (ary instanceof ASTId) ary = Env.staticLookup((ASTId)ary);
    // Get the trim
    try {
      _trim = ((ASTNum) (E.skipWS().parse())).dbl();
    } catch (ClassCastException e) {
      e.printStackTrace();
      throw new IllegalArgumentException("Argument `trim` expected to be a number.");
    }
    // Get the na.rm
    AST a = E._env.lookup((ASTId)E.skipWS().parse());
    try {
      _narm = ((ASTNum) a).dbl() == 1;
    } catch (ClassCastException e) {
      e.printStackTrace();
      throw new IllegalArgumentException("Argument `na.rm` expected to be a number.");
    }
    // Finish the rest
    ASTMean res = (ASTMean) clone();
    res._asts = new AST[]{ary};
    return res;
  }

  @Override void exec(Env e, AST arg1, AST[] args) {
    arg1.exec(e);
    e._global._frames.put(Key.make().toString(), e.peekAry());
    if (args != null) {
      if (args.length > 2) throw new IllegalArgumentException("Too many arguments passed to `mean`");
      for (AST a : args) {
        if (a instanceof ASTId) {
          _narm = ((ASTNum) e.lookup((ASTId) a)).dbl() == 1;
        } else if (a instanceof ASTNum) {
          _trim = ((ASTNum) a).dbl();
        }
      }
    }
    apply(e);
  }

  @Override protected void apply(Env env) {
    if (env.isNum()) return;
    Frame fr = env.popAry(); // get the frame w/o sub-reffing
    if (fr.numCols() > 1 && fr.numRows() > 1)
      throw new IllegalArgumentException("mean does not apply to multiple cols.");
    for (Vec v : fr.vecs()) if (v.isEnum())
      throw new IllegalArgumentException("mean only applies to numeric vector.");
    if (fr.numCols() > 1) {
      double mean=0;
      for(Vec v : fr.vecs()) mean += v.at(0);
      env.push(new ValNum(mean/fr.numCols()));
    } else {
      MeanNARMTask t = new MeanNARMTask(_narm).doAll(fr.anyVec()).getResult();
      if (t._rowcnt == 0 || Double.isNaN(t._sum)) {
        double ave = Double.NaN;
        env.push(new ValNum(ave));
      } else {
        double ave = t._sum / t._rowcnt;
        env.push(new ValNum(ave));
      }
    }
  }

  @Override double[] map(Env e, double[] in, double[] out, AST[] args) {
    if (args != null) {
      if (args.length > 2) throw new IllegalArgumentException("Too many arguments passed to `mean`");
      for (AST a : args) {
        if (a instanceof ASTId) {
          _narm = ((ASTNum) e.lookup((ASTId) a)).dbl() == 1;
        } else if (a instanceof ASTNum) {
          _trim = ((ASTNum) a).dbl();
        }
      }
    }
    if (out == null || out.length < 1) out = new double[1];
    double s = 0;  int cnt=0;
    for (double v : in) if( !Double.isNaN(v) ) { s+=v; cnt++; }
    out[0] = s/cnt;
    return out;
  }

  static class MeanNARMTask extends MRTask<MeanNARMTask> {
    // IN
    boolean _narm;    // remove NAs
    double  _trim;    // trim each end of the column -- unimplemented: requires column sort
    int     _nrow;    // number of rows in the colun -- useful for trim

    // OUT
    long   _rowcnt;
    double _sum;
   MeanNARMTask(boolean narm) {
     _narm = narm;
//     _trim = trim;
//     _nrow = nrow;
//     if (_trim != 0) {
//       _start = (long)Math.floor(_trim * (nrow - 1));
//       _end   = (long)(nrow - Math.ceil(_trim * (nrow - 1)));
//     }
   }
    @Override public void map(Chunk c) {
      if (c.vec().isEnum() || c.vec().isUUID()) { _sum = Double.NaN; _rowcnt = 0; return;}
      if (_narm) {
        for (int r = 0; r < c._len; r++)
          if (!c.isNA(r)) { _sum += c.atd(r); _rowcnt++;}
      } else {
        for (int r = 0; r < c._len; r++)
          if (c.isNA(r)) { _rowcnt = 0; _sum = Double.NaN; return; } else { _sum += c.atd(r); _rowcnt++; }
      }
    }
    @Override public void reduce(MeanNARMTask t) {
      _rowcnt += t._rowcnt;
      _sum += t._sum;
    }
  }
}

class ASTTable extends ASTUniPrefixOp {
  ASTTable() { super(new String[]{"table", "..."}); }
  @Override protected String opStr() { return "table"; }
  @Override protected ASTOp make() { return new ASTTable(); }

  @Override
  protected ASTTable parse_impl(Exec E) {
    AST ary = E.parse();
    if (ary instanceof ASTId) ary = Env.staticLookup((ASTId)ary);
    AST two = E.skipWS().parse();
    if (two instanceof ASTString) two = new ASTNull();
    if (two instanceof ASTId) two = Env.staticLookup((ASTId)two);
    ASTTable res = (ASTTable)clone();
    res._asts = new AST[]{ary, two}; //two is pushed on, then ary is pushed on
    return res;
  }

  @Override protected void apply(Env env) {
    Frame two = env.peekType() == Env.NULL ? null : env.popAry();
    if (two == null) env.pop();
    Frame one = env.popAry();

    // Rules: two != null => two.numCols == one.numCols == 1
    //        two == null => one.numCols == 1 || one.numCols == 2
    // Anything else is IAE

    if (two != null)
      if (two.numCols() != 1 || one.numCols() != 1)
        throw new IllegalArgumentException("`table` supports at *most* two vectors");
      else if (one.numCols() < 1 || one.numCols() > 2)
        throw new IllegalArgumentException("`table` supports at *most* two vectors and at least one vector.");

    Frame fr;
    if (two != null) fr = new Frame(one.add(two));
    else fr = one;

    final int ncol;
    if ((ncol = fr.vecs().length) > 2)
      throw new IllegalArgumentException("table does not apply to more than two cols.");
    for (int i = 0; i < ncol; i++)
      if (!fr.vecs()[i].isInt())
        throw new IllegalArgumentException("table only applies to integer vectors.");
    Vec dataLayoutVec;
    Frame fr2;
    String colnames[];
    String[][] d = new String[ncol+1][];

    if (ncol == 1) {
      final int min = (int)fr.anyVec().min();
      final int max = (int)fr.anyVec().max();
      colnames = new String[]{fr.name(0), "Count"};
      d[0] = fr.anyVec().domain(); // should always be null for all neg values!
      d[1] = null;

      // all pos
      if (min >= 0) {
        UniqueColumnCountTask t = new UniqueColumnCountTask(max,false,false,0).doAll(fr.anyVec());
        final long[] cts = t._cts;
        dataLayoutVec = Vec.makeCon(0, cts.length);

        // second pass to build the result frame
        fr2 = new MRTask() {
          @Override public void map(Chunk[] c, NewChunk[] cs) {
            for (int i = 0; i < c[0]._len; ++i) {
              int idx = (int) (i + c[0].start());
              if (cts[idx] == 0) continue;
              cs[0].addNum(idx);
              cs[1].addNum(cts[idx]);
            }
          }
        }.doAll(2, dataLayoutVec).outputFrame(colnames, d);

        // all neg  -- flip the sign and count...
      } else if (min <= 0 && max <= 0) {
        UniqueColumnCountTask t = new UniqueColumnCountTask(-1*min,true,false,0).doAll(fr.anyVec());
        final long[] cts = t._cts;
        dataLayoutVec = Vec.makeCon(0, cts.length);

        // second pass to build the result frame
        fr2 = new MRTask() {
          @Override public void map(Chunk[] c, NewChunk[] cs) {
            for (int i = 0; i < c[0]._len; ++i) {
              int idx = (int) (i + c[0].start());
              if (cts[idx] == 0) continue;
              cs[0].addNum(idx * -1);
              cs[1].addNum(cts[idx]);
            }
          }
        }.doAll(2, dataLayoutVec).outputFrame(colnames, d);

        // mixed
      } else {
        UniqueColumnCountTask t = new UniqueColumnCountTask(max+-1*min,false,true,max).doAll(fr.anyVec()); // pivot around max value... vals > max are negative
        final long[] cts = t._cts;
        dataLayoutVec = Vec.makeCon(0, cts.length);

        // second pass to build the result frame
        fr2 = new MRTask() {
          @Override public void map(Chunk[] c, NewChunk[] cs) {
            for (int i = 0; i < c[0]._len; ++i) {
              int idx = (int) (i + c[0].start());
              if (cts[idx] == 0) continue;
              cs[0].addNum(idx > max ? (idx-max)*-1 : idx);
              cs[1].addNum(cts[idx]);
            }
          }
        }.doAll(2, dataLayoutVec).outputFrame(colnames, d);
      }

    } else {

      // 2 COLUMN CASE
      colnames = new String[]{fr.name(0), fr.name(1), "count"};
      long s = System.currentTimeMillis();

      // Strategy: Avoid doing NBHM reduces (which is way too slow).
      //  1. Build NBHS of Groups (useful for
      Uniq2ColTsk u = new Uniq2ColTsk().doAll(fr);
      Log.info("Finished gathering uniq groups in: " + (System.currentTimeMillis() - s) / 1000. + " (s)");

      final long[] pairs = new long[u._s.size()];
      int i=0;
      for (Object o : u._s) pairs[i++] = (long) o;
      dataLayoutVec = Vec.makeCon(0, pairs.length);

      s = System.currentTimeMillis();
      NewHashMap h = new NewHashMap(pairs).doAll(dataLayoutVec);
      Log.info("Finished creating new HashMap in: " + (System.currentTimeMillis() - s) / 1000. + " (s)");

      s = System.currentTimeMillis();
      final long[] cnts = new CountUniq2ColTsk(h._s).doAll(fr)._cnts;
      Log.info("Finished gathering counts in: " + (System.currentTimeMillis() - s) / 1000. + " (s)");

      d[0] = fr.vec(0).domain();
      d[1] = fr.vec(1).domain();

      fr2 = new MRTask() {
        @Override
        public void map(Chunk[] c, NewChunk[] cs) {
          int start = (int)c[0].start();
          for (int i = 0; i < c[0]._len; ++i) {
            long[] g = unmix(pairs[i+start]);
            cs[0].addNum(g[0]);
            cs[1].addNum(g[1]);
            cs[2].addNum(cnts[i+start]);
          }
        }
      }.doAll(ncol + 1, dataLayoutVec).outputFrame(colnames, d);
    }
    Keyed.remove(dataLayoutVec._key);
    env.pushAry(fr2);
  }

  // gets vast majority of cases and is stupidly fast (35x faster than using UniqueTwoColumnTask)
  private static class UniqueColumnCountTask extends MRTask<UniqueColumnCountTask> {
    long[] _cts;
    final int _max;
    final boolean _flip;
    final boolean _mixed;
    final int _piv;
    public UniqueColumnCountTask(int max, boolean flip, boolean mixed, int piv) { _max = max; _flip = flip; _mixed = mixed; _piv = piv; }
    @Override public void map( Chunk c ) {
      _cts = MemoryManager.malloc8(_max+1);
      // choose the right hot loop
      if (_flip) {
        for (int i = 0; i < c._len; ++i) {
          if (c.isNA(i)) continue;
          int val = (int) (-1 * c.at8(i));
          _cts[val]++;
        }
      } else if (_mixed) {
        for (int i = 0; i < c._len; ++i) {
          if (c.isNA(i)) continue;
          int val = (int) (c.at8(i));
          int idx = val < 0 ? -1*val + _piv : val;
          _cts[idx]++;
        }
      } else {
        for (int i = 0; i < c._len; ++i) {
          if (c.isNA(i)) continue;
          int val = (int) (c.at8(i));
          _cts[val]++;
        }
      }
    }
    @Override public void reduce(UniqueColumnCountTask t) { ArrayUtils.add(_cts, t._cts); }
  }

  private static class Uniq2ColTsk extends MRTask<Uniq2ColTsk> {
    NonBlockingHashSet<Long> _s;
    @Override public void setupLocal() { _s = new NonBlockingHashSet<>(); }
    @Override public void map(Chunk[] c) {
      for (int i=0;i<c[0]._len;++i)
        _s.add(mix(c[0].at8(i), c[1].at8(i)));
    }
    @Override public void reduce(Uniq2ColTsk t) { if (_s!=t._s) _s.addAll(t._s); }

    @Override public AutoBuffer write_impl( AutoBuffer ab ) {
      if( _s == null ) return ab.put4(0);
      ab.put4(_s.size());
      for( Long g : _s ) {ab.put8(g); }
      return ab;
    }

    @Override public Uniq2ColTsk read_impl( AutoBuffer ab ) {
      int len = ab.get4();
      if( len == 0 ) return this;
      _s = new NonBlockingHashSet<>();
      for( int i=0; i<len; i++ ) { _s.add(ab.get8());}
      return this;
    }
  }

  /** http://szudzik.com/ElegantPairing.pdf */
  private static long mix(long A, long B) {
    long a=A,b=B;
    long a1 = (a<<=1) >= 0 ? a : -1 * a - 1;
    long b1 = (b<<=1) >= 0 ? b : -1 * b - 1;
    long v = (a1 >= b1 ? a1 * a1 + a1 + b1 : a1 + b1 * b1) >> 1; // pairing fcn
    return a < 0 && b < 0 || a >= 0 && b >= 0 ? v : -v - 1;
  }

  // always returns a long[] of length 2;
  // long[0] -> A, long[1] -> B
  private static long[] unmix(long z) {
    long rflr = (long)Math.floor(Math.sqrt(z));
    long rflr2=rflr*rflr;
    long z_rflr2 = z-rflr2;
    return z_rflr2 < rflr ?  new long[]{z_rflr2,rflr} : new long[]{rflr, z_rflr2-rflr};
  }

  private static class NewHashMap extends MRTask<NewHashMap> {
    NonBlockingHashMap<Long, Integer> _s;
    final long[] _m;
    NewHashMap(long[] m) { _m = m; }
    @Override public void setupLocal() { _s = new NonBlockingHashMap<>();}
    @Override public void map(Chunk[] c) {
      int start = (int)c[0].start();
      for (int i = 0; i < c[0]._len; ++i)
        _s.put(_m[i + start], i+start);
    }
    @Override public void reduce(NewHashMap t) { if (_s != t._s) _s.putAll(t._s); }

    @Override public AutoBuffer write_impl( AutoBuffer ab ) {
      if( _s == null ) return ab.put4(0);
      ab.put4(_s.size());
      for( Long l : _s.keySet() ) {ab.put8(l); ab.put4(_s.get(l)); }
      return ab;
    }

    @Override public NewHashMap read_impl( AutoBuffer ab ) {
      int len = ab.get4();
      if( len == 0 ) return this;
      _s = new NonBlockingHashMap<>();
      for( int i=0;i<len;i++ ) _s.put(ab.get8(), ab.get4());
      return this;
    }
  }

  private static class CountUniq2ColTsk extends MRTask<CountUniq2ColTsk> {
    private static final Unsafe _unsafe = UtilUnsafe.getUnsafe();
    final NonBlockingHashMap<Long, Integer> _m;
    // out
    long[] _cnts;
    private static final int _b = _unsafe.arrayBaseOffset(long[].class);
    private static final int _s = _unsafe.arrayIndexScale(long[].class);
    private static long ssid(int i) { return _b + _s*i; } // Scale and Shift

    CountUniq2ColTsk(NonBlockingHashMap<Long, Integer> s) {_m = s; }
    @Override public void setupLocal() { _cnts = MemoryManager.malloc8(_m.size()); }
    @Override public void map(Chunk[] cs) {
      for (int i=0; i < cs[0]._len; ++i) {
        int h = _m.get(mix(cs[0].at8(i), cs[1].at8(i)));
        long offset = ssid(h);
        long c = _cnts[h];
        while(!_unsafe.compareAndSwapLong(_cnts,offset,c,c+1))  //yee-haw
          c = _cnts[h];
      }
    }
    @Override public void reduce(CountUniq2ColTsk t) { if (_cnts != t._cnts) ArrayUtils.add(_cnts, t._cnts); }
  }
}

// Conditional merge of two Frames/Vecs
// Result is always the height as "tst"
// That means we support R semantics here and do the replication of true/false as needed.
// We also do the same thing R does with factor levels -> replace with their int value...
// What we do NOT support is the case when true.numCols() != false.numCols()
// The pseudo-code is: [t <- true; t[!tst] < false[!tst]; t]

// Cases to consider: (2^2 possible input types f: Frame, v: Vec (1D frame), the tst may only be a single column bit vec.
//
//       tst | true | false
//       ----|------|------
//     1.  v     f      f
//     2.  v     f      v
//     3.  v     v      f
//     4.  v     v      v
//
//  Additionally, how to cut/expand frames `true` and `false` to match `tst`.

class ASTIfElse extends ASTUniPrefixOp {
  static final String VARS[] = new String[]{"ifelse","tst","true","false"};

  ASTIfElse( ) { super(VARS); }
  @Override protected ASTOp make() {return new ASTIfElse();}
  @Override protected String opStr() { return "ifelse"; }
  @Override
  protected ASTIfElse parse_impl(Exec E) {
    AST tst = E.parse();
    if (tst instanceof ASTId) tst = Env.staticLookup((ASTId) tst);

    // still have an instance of ASTId, and lookup gives 0 (%FALSE) or 1 (%TRUE)
    if (tst instanceof ASTId) {
      try {
        double d = ((ASTNum) E._env.lookup((ASTId)tst))._d;
        if (d == 0 || d == 1) {  // FALSE or TRUE
          tst = new ASTFrame(new Frame(Key.make(), null, new Vec[]{Vec.makeCon(d, 1)} ) );
        }
      } catch (ClassCastException e) {
        throw new IllegalArgumentException("`test` must be a frame or TRUE/FALSE");
      }
    }
    AST yes = E.skipWS().parse(); // could be num
    if (yes instanceof ASTId) yes = Env.staticLookup((ASTId)yes);
    AST no  = E.skipWS().parse(); // could be num
    if (no instanceof ASTId) no = Env.staticLookup((ASTId)no);
    ASTIfElse res = (ASTIfElse)clone();
    res._asts = new AST[]{no,yes,tst};
    return res;
  }

  // return frame compatible to tgt
  private Frame adaptToTst(Frame src, Frame tgt) {
    Key k = src._key == null ? Key.make() : src._key;
    // need to pute src in DKV if not in there
    if (src._key == null || DKV.get(src._key) == null)
      DKV.put(k, new Frame(k,src.names(),src.vecs()));

    // extend src
    StringBuilder sb=null;
    if (src.numRows() < tgt.numRows()) {
      // rbind the needed rows
      int nrbins = 1 + (int)((tgt.numRows() - src.numRows()) / src.numRows());
      long remainder = tgt.numRows() % src.numRows();
      sb = new StringBuilder("(rbind ");
      for (int i = 0; i < nrbins; ++i) sb.append("%").append(k).append((i == (nrbins - 1) && remainder<0) ? "" : " ");
      sb.append(remainder > 0 ? "([ %"+k+" (: #0 #"+(remainder-1)+") \"null\"))" : ")");
      Log.info("extending frame:" + sb.toString());

    // reduce src
    } else if (src.numRows() > tgt.numRows()) {
      long rmax = tgt.numRows() - 1;
      sb = new StringBuilder("([ %"+k+" (: #0 #"+rmax+"))");
    }

    if (sb != null) {
      Env env=null;
      Frame res;
      try {
        env = Exec.exec(sb.toString());
        res = env.popAry();
        res.unlock_all();
      } catch (Exception e) {
        throw new H2OIllegalArgumentException("Bad expression Rapids: " + sb.toString(), "Bad expression Rapids: " + sb.toString() + "; exception: " + e.toString());
      } finally {
        if (env!=null)env.unlock();
      }
      Frame ret = tgt.makeCompatible(res);
//      if (env != null) env.cleanup(ret==res?null:res, (Frame)DKV.remove(k).get());
      return ret;
    }
    src = DKV.remove(k).get();
    Frame ret = tgt.makeCompatible(src);
    if (src != ret) src.delete();
    return ret;
  }

  private Frame adaptToTst(double d, Frame tgt) {
    Frame v = new Frame(Vec.makeCon(d, tgt.numRows()));
    Frame ret = tgt.makeCompatible(v);
    if (ret != v) v.delete();
    return ret;
  }

  @Override protected void apply(Env env) {
    if (!env.isAry()) throw new IllegalArgumentException("`test` argument must be a frame: ifelse(`test`, `yes`, `no`)");
    Frame tst = env.popAry();
    if (tst.numCols() != 1)
      throw new IllegalArgumentException("`test` has "+tst.numCols()+" columns. `test` must have exactly 1 column.");
    Frame yes=null; double dyes=0;
    Frame no=null; double dno=0;
    if (env.isAry()) yes = env.popAry(); else dyes = env.popDbl();
    if (env.isAry()) no  = env.popAry(); else dno  = env.popDbl();

    if (yes != null && no != null) {
      if (yes.numCols() != no.numCols())
        throw new IllegalArgumentException("Column mismatch between `yes` and `no`. `yes` has" + yes.numCols() + "; `no` has " + no.numCols() + ".");
    } else if (yes != null) {
      if (yes.numCols() != 1)
        throw new IllegalArgumentException("Column mismatch between `yes` and `no`. `yes` has" + yes.numCols() + "; `no` has " + 1 + ".");
    } else if (no != null) {
      if (no.numCols() != 1)
        throw new IllegalArgumentException("Column mismatch between `yes` and `no`. `yes` has" + 1 + "; `no` has " + no.numCols() + ".");
    }

    Frame a_yes = yes == null ? adaptToTst(dyes, tst) : adaptToTst(yes,tst);
    Frame a_no  = no == null ? adaptToTst(dno, tst) : adaptToTst(no, tst);
    Frame frtst = (new Frame(tst)).add(a_yes).add(a_no);
    final int ycols = a_yes.numCols();

    // Run a selection picking true/false across the frame
    Frame fr2 = new MRTask() {
      @Override public void map( Chunk chks[], NewChunk nchks[] ) {
        int rows = chks[0]._len;
        int cols = chks.length;
        Chunk pred = chks[0];
        for (int r=0;r < rows;++r) {
          for (int c = (pred.atd(r) != 0 ? 1 : ycols + 1), col = 0; c < (pred.atd(r) != 0 ?ycols+1:cols); ++c) {
            if (chks[c].vec().isUUID())
              nchks[col++].addUUID(chks[c], r);
            else if (chks[c].vec().isString())
              nchks[col++].addStr(chks[c].atStr(new ValueString(), r));
            else
              nchks[col++].addNum(chks[c].atd(r));
          }
        }
      }
    }.doAll(yes==null?1:yes.numCols(),frtst).outputFrame(yes==null?(new String[]{"C1"}):yes.names(),null/*same as R: no domains*/);
    env.pushAry(fr2);
  }
}

class ASTCut extends ASTUniPrefixOp {
  String[] _labels = null;
  double[] _cuts;
  boolean _includelowest = false;
  boolean _right = true;
  double _diglab = 3;
  ASTCut() { super(new String[]{"cut", "ary", "breaks", "labels", "include.lowest", "right", "dig.lab"});}
  @Override protected String opStr() { return "cut"; }
  @Override protected ASTOp make() {return new ASTCut();}
  protected ASTCut parse_impl(Exec E) {
    AST ary = E.parse();
    if (ary instanceof ASTId) ary = Env.staticLookup((ASTId)ary);
    // breaks first
    String[] cuts;
    try {
      cuts = E.skipWS().peek() == '{'
              ? E.xpeek('{').parseString('}').split(";")
              : E.peek() == '#' ? new String[]{Double.toString(((ASTNum) E.parse()).dbl())}
              : new String[]{E.parseString(E.peekPlus())};
    } catch (ClassCastException e) {
      e.printStackTrace();
      throw new IllegalArgumentException("Argument `breaks` was malformed. Bad AST input.");
    }
    for (int i = 0; i < cuts.length; ++i) cuts[i] = cuts[i].replace("\"", "").replace("\'", "");
    _cuts = new double[cuts.length];
    for (int i = 0; i < cuts.length; ++i) _cuts[i] = Double.valueOf(cuts[i]);
    // labels second
    try {
      _labels = E.skipWS().peek() == '{' ? E.xpeek('{').parseString('}').split(";") : new String[]{E.parseString(E.peekPlus())};
    } catch (ClassCastException e) {
      e.printStackTrace();
      throw new IllegalArgumentException("Argument `labels` was malformed. Bad AST input.");
    }
    // cleanup _labels
    for (int i = 0; i < _labels.length; ++i) _labels[i] = _labels[i].replace("\"", "").replace("\'", "");
    if (_labels.length==1 && _labels[0].equals("null")) _labels = null;
    AST inc_lowest = E.skipWS().parse();
    inc_lowest = E._env.lookup((ASTId)inc_lowest);
    try {
      _includelowest = ((ASTNum) inc_lowest).dbl() == 1;
    } catch (ClassCastException e) {
      e.printStackTrace();
      throw new IllegalArgumentException("Argument `include.lowest` expected to be TRUE/FALSE.");
    }
    AST right = E.skipWS().parse();
    right = E._env.lookup((ASTId)right);
    try {
      _right = ((ASTNum) right).dbl() == 1;
    } catch (ClassCastException e) {
      e.printStackTrace();
      throw new IllegalArgumentException("Argument `right` expected to be a TRUE/FALSE.");
    }
    ASTNum diglab;
    try {
      diglab = (ASTNum) E.skipWS().parse();
    } catch (ClassCastException e) {
      e.printStackTrace();
      throw new IllegalArgumentException("Argument `dig.lab` expected to be a number.");
    }
    _diglab = diglab.dbl();
    _diglab = _diglab >= 12 ? 12 : _diglab; // cap at 12 digits
    ASTCut res = (ASTCut) clone();
    res._asts = new AST[]{ary};
    return res;
  }

  private String left() { return _right ? "(" : "["; }
  private String rite() { return _right ? "]" : ")"; }
  @Override protected void apply(Env env) {
    Frame fr = env.popAry();
    if(fr.vecs().length != 1 || fr.vecs()[0].isEnum())
      throw new IllegalArgumentException("First argument must be a numeric column vector");

    double fmin = fr.anyVec().min();
    double fmax = fr.anyVec().max();

    int nbins = _cuts.length - 1;  // c(0,10,100) -> 2 bins (0,10] U (10, 100]
    double width;
    if (nbins == 0) {
      if (_cuts[0] < 2) throw new IllegalArgumentException("The number of cuts must be >= 2. Got: "+_cuts[0]);
      // in this case, cut the vec into _cuts[0] many pieces of equal length
      nbins = (int) Math.floor(_cuts[0]);
      width = (fmax - fmin)/nbins;
      _cuts = new double[nbins];
      _cuts[0] = fmin - 0.001*(fmax - fmin);
      for (int i = 1; i < _cuts.length; ++i) _cuts[i] = (i == _cuts.length-1) ? (fmax + 0.001*(fmax-fmin))  : (fmin + i*width);
    }
    width = (fmax - fmin)/nbins;
    if(width == 0) throw new IllegalArgumentException("Data vector is constant!");
    if (_labels != null && _labels.length != nbins) throw new IllegalArgumentException("`labels` vector does not match the number of cuts.");

    // Construct domain names from _labels or bin intervals if _labels is null
    final double cuts[] = _cuts;

    // first round _cuts to dig.lab decimals: example floor(2.676*100 + 0.5) / 100
    for (int i = 0; i < _cuts.length; ++i) _cuts[i] = Math.floor(_cuts[i] * Math.pow(10,_diglab) + 0.5) / Math.pow(10,_diglab);

    String[][] domains = new String[1][nbins];
    if (_labels == null) {
      domains[0][0] = (_includelowest ? "[" : left()) + _cuts[0] + "," + _cuts[1] + rite();
      for (int i = 1; i < (_cuts.length - 1); ++i)  domains[0][i] = left() + _cuts[i] + "," + _cuts[i+1] + rite();
    } else domains[0] = _labels;

    final boolean incLow = _includelowest;
    Frame fr2 = new MRTask() {
      @Override public void map(Chunk c, NewChunk nc) {
        int rows = c._len;
        for (int r = 0; r < rows; ++r) {
          double x = c.atd(r);
          if (Double.isNaN(x) || (incLow  && x <  cuts[0])
                              || (!incLow && x <= cuts[0])
                              || (_right  && x >  cuts[cuts.length-1])
                              || (!_right && x >= cuts[cuts.length-1])) nc.addNum(Double.NaN);
          else {
            for (int i = 1; i < cuts.length; ++i) {
              if (_right) {
                if (x <= cuts[i]) {
                  nc.addNum(i - 1);
                  break;
                }
              } else if (x < cuts[i]) { nc.addNum(i-1); break; }
            }
          }
        }
      }
    }.doAll(1, fr).outputFrame(fr.names(), domains);
    env.pushAry(fr2);
  }
}

class ASTAsNumeric extends ASTUniPrefixOp {
  ASTAsNumeric() { super(new String[]{"as.numeric", "ary"}); }
  @Override protected String opStr() { return "as.numeric"; }
  @Override protected ASTOp make() {return new ASTAsNumeric(); }
  protected ASTAsNumeric parse_impl(Exec E) {
    AST ary = E.parse();
    if (ary instanceof ASTId) ary = Env.staticLookup((ASTId)ary);
    ASTAsNumeric res = (ASTAsNumeric) clone();
    res._asts = new AST[]{ary};
    return res;
  }
  @Override protected void apply(Env env) {
    Frame ary = env.peekAry();
    Vec[] nvecs = new Vec[ary.numCols()];
    for (int c = 0; c < ary.numCols(); ++c)
      nvecs[c] = ary.vecs()[c].toInt();
    Frame v = new Frame(ary._names, nvecs);
    env.poppush(1, new ValFrame(v));
  }
}

class ASTFactor extends ASTUniPrefixOp {
  ASTFactor() { super(new String[]{"", "ary"});}
  @Override protected String opStr() { return "as.factor"; }
  @Override protected ASTOp make() {return new ASTFactor();}
  protected ASTFactor parse_impl(Exec E) {
    AST ary = E.parse();
    if (ary instanceof ASTId) ary = Env.staticLookup((ASTId)ary);
    ASTFactor res = (ASTFactor) clone();
    res._asts = new AST[]{ary};
    return res;
  }
  @Override protected void apply(Env env) {
    Frame ary = env.popAry();
    if( ary.numCols() != 1 ) throw new IllegalArgumentException("factor requires a single column");
    Vec v0 = ary.anyVec();
    if( v0.isEnum() ) {
      env.pushAry(ary);
      return;
    }
    Vec v1 = v0.toEnum(); // toEnum() creates a new vec --> must be cleaned up!
    Frame fr = new Frame(ary._names, new Vec[]{v1});
    env.pushAry(fr);
  }
}

class ASTCharacter extends ASTUniPrefixOp {
  ASTCharacter() { super(new String[]{"", "ary"});}
  @Override protected String opStr() { return "as.character"; }
  @Override protected ASTOp make() {return new ASTFactor();}
  protected ASTCharacter parse_impl(Exec E) {
    AST ary = E.parse();
    if (ary instanceof ASTId) ary = Env.staticLookup((ASTId)ary);
    ASTCharacter res = (ASTCharacter) clone();
    res._asts = new AST[]{ary};
    return res;
  }
  @Override protected void apply(Env env) {
    Frame ary = env.popAry();
    if( ary.numCols() != 1 ) throw new IllegalArgumentException("character requires a single column");
    Vec v0 = ary.anyVec();
    Vec v1 = v0.isString() ? null : v0.toStringVec(); // toEnum() creates a new vec --> must be cleaned up!
    Frame fr = new Frame(ary._names, new Vec[]{v1 == null ? v0.makeCopy(null) : v1});
    env.pushAry(fr);
  }
}

// Variable length; flatten all the component arys
class ASTCat extends ASTUniPrefixOp {
  // form of this is (c ASTSpan)
  @Override protected String opStr() { return "c"; }
  public ASTCat( ) { super(new String[]{"cat","dbls", "..."});}
  @Override protected ASTOp make() {return new ASTCat();}
  @Override
  protected ASTCat parse_impl(Exec E) {
    ASTSeries a;
    try {
      if (!E.hasNext()) throw new IllegalArgumentException("End of input unexpected. Badly formed AST.");
      a = (ASTSeries) E.parse();
    } catch (ClassCastException e) {
      throw new IllegalArgumentException("Expected ASTSeries object. Badly formed AST.");
    }

    ASTCat res = (ASTCat) clone();
    res._asts = new AST[]{a};
    return res;
  }

  @Override protected void apply(Env env) {
    final ValSeries s = (ValSeries) env.pop();
    int id_span =0;
    long len = s._idxs.length;
    if (s._spans != null) {
      for (ASTSpan as : s._spans) len += (as._max - as._min + 1);
    }
    // now make an mapping of ValSeries -> Vec indices
    ArrayList<Long> idxs = new ArrayList<>();
    ArrayList<ASTSpan> spans = new ArrayList<>();
    long cur_id=0;
    for (int o : s._order) {
      if (o == 0) { // span
        assert s._spans != null;
        long id_min = cur_id;
        long id_max = cur_id + s._spans[id_span]._max - s._spans[id_span]._min;
        cur_id+=(s._spans[id_span]._max-s._spans[id_span++]._min+1);
        spans.add(new ASTSpan(id_min, id_max));
      } else {      // idx
        idxs.add(cur_id++);
      }
    }
    long[] idxsl = new long[idxs.size()];
    for (int i =0; i < idxsl.length; ++ i) idxsl[i] = idxs.get(i);
    final ValSeries ids = new ValSeries(idxsl, spans.toArray(new ASTSpan[spans.size()]));
    ids._order = s._order;

    Frame fr = new MRTask() {
      @Override public void map(Chunk[] cs) {
        Chunk c = cs[0];
        for (int r = 0; r < c._len; ++r) {
          long cur = c.start() + r;
          c.set(r, maprow(cur, ids, s));
        }
      }
    }.doAll(Vec.makeZero(len))._fr;

    env.pushAry(fr);
  }

  private long maprow(long cur, ValSeries ids, ValSeries s) {
    // get the location of the id in ids. This maps to the value in s.
    int span_idx = -1;
    int idxs_idx = 0;
    long at_value = -1;
    for (int o : ids._order) {
      if (o == 0) {  // span
        if (ids._spans[++span_idx].contains(cur)) {
          at_value = cur - ids._spans[span_idx]._min;
          break;
        }
      } else {
        boolean _br = false;
        if (idxs_idx >= 0) {
          for (int i = 0; i < ids._idxs.length; ++i) {
            if (ids._idxs[i] == cur) {
              at_value = i;
              _br = true;
              break;
            }
          }
          if (_br) break; else --idxs_idx;
        }
      }
    }
    if (span_idx >= 0) {
      return s._spans[span_idx]._min + at_value;
    } else return s._idxs[(int)at_value];
  }
}

//class ASTFindInterval extends ASTUniPrefixOp {
//  protected static boolean _rclosed;
//  protected static double _x;
//
//  ASTFindInterval() { super(new String[]{"findInterval", "x", "vec", "rightmost.closed"}); }
//  @Override protected String opStr() { return "findInterval"; }
//  @Override protected ASTOp make() { return new ASTFindInterval(); }
//  @Override ASTFindInterval parse_impl(Exec E) {
//    // First argument must be a num, anything else throw IAE
//    AST x = E.skipWS().parse();
//    if (! (x instanceof ASTNum)) throw new IllegalArgumentException("First argument to findInterval must be a single number. Got: " + x.toString());
//    _x =  ((ASTNum)(E.skipWS().parse())).dbl();
//    // Get the ary
//    AST ary = E.parse();
//    // Get the rightmost.closed
//    AST a = E._env.lookup((ASTId)E.skipWS().parse());
//    _rclosed = ((ASTNum)a).dbl() == 1;
//    // Finish the rest
//    ASTFindInterval res = (ASTFindInterval) clone();
//    res._asts = new AST[]{ary};
//    return res;
//  }
//
//  @Override protected void apply(Env env) {
//    final boolean rclosed = _rclosed;
//
//    if(env.isNum()) {
//      final double cutoff = _x;
//
//      Frame fr = env.popAry();
//      if(fr.numCols() != 1 || fr.vecs()[0].isEnum())
//        throw new IllegalArgumentException("Argument must be a single numeric column vector. Got an array with " + fr.numCols() + " columns. Column was an enum: " + fr.vecs()[0].isEnum());
//
//      Frame fr2 = new MRTask() {
//        @Override public void map(Chunk chk, NewChunk nchk) {
//          for(int r = 0; r < chk._len; r++) {
//            double x = chk.atd(r);
//            if(Double.isNaN(x))
//              nchk.addNum(Double.NaN);
//            else {
//              if(rclosed)
//                nchk.addNum(x > cutoff ? 1 : 0);   // For rightmost.closed = TRUE
//              else
//                nchk.addNum(x >= cutoff ? 1 : 0);
//            }
//          }
//        }
//      }.doAll(1,fr).outputFrame(fr._names, fr.domains());
//      env.subRef(fr, skey);
//      env.pop();
//      env.push(fr2);
//    } else if(env.isAry()) {
//      Frame ary = env.popAry();
//      String skey1 = env.key();
//      if(ary.vecs().length != 1 || ary.vecs()[0].isEnum())
//        throw new IllegalArgumentException("Second argument must be a numeric column vector");
//      Vec brks = ary.vecs()[0];
//      // TODO: Check that num rows below some cutoff, else this will likely crash
//
//      // Check if vector of cutoffs is sorted in weakly ascending order
//      final int len = (int)brks.length();
//      final double[] cutoffs = new double[len];
//      for(int i = 0; i < len-1; i++) {
//        if(brks.at(i) > brks.at(i+1))
//          throw new IllegalArgumentException("Second argument must be sorted in non-decreasing order");
//        cutoffs[i] = brks.at(i);
//      }
//      cutoffs[len-1] = brks.at(len-1);
//
//      Frame fr = env.popAry();
//      String skey2 = env.key();
//      if(fr.vecs().length != 1 || fr.vecs()[0].isEnum())
//        throw new IllegalArgumentException("First argument must be a numeric column vector");
//
//      Frame fr2 = new MRTask() {
//        @Override public void map(Chunk chk, NewChunk nchk) {
//          for(int r = 0; r < chk._len; r++) {
//            double x = chk.atd(r);
//            if(Double.isNaN(x))
//              nchk.addNum(Double.NaN);
//            else {
//              double n = Arrays.binarySearch(cutoffs, x);
//              if(n < 0) nchk.addNum(-n-1);
//              else if(rclosed && n == len-1) nchk.addNum(n);   // For rightmost.closed = TRUE
//              else nchk.addNum(n+1);
//            }
//          }
//        }
//      }.doAll(1,fr).outputFrame(fr._names, fr.domains());
//      env.subRef(ary, skey1);
//      env.subRef(fr, skey2);
//      env.pop();
//      env.push(fr2);
//    }
//  }
//}

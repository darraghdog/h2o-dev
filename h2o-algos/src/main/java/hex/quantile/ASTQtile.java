package hex.quantile;

import water.H2O;
import water.rapids.ASTUniPrefixOp;
import water.rapids.Exec;
import water.rapids.Env;

// Compute exact quantiles given a set of cutoffs, using multipass binning algo.
class ASTQtile extends ASTUniPrefixOp {
  protected static double[] _probs = null;  // if probs is null, pop the _probs frame etc.

  @Override protected String opStr() { return "quantile"; }

  public ASTQtile() { super(new String[]{"quantile","x","probs"});}
  @Override protected ASTQtile make() { return new ASTQtile(); }
  @Override protected ASTQtile parse_impl(Exec E) {
    throw H2O.unimpl();
    //// Get the ary
    //AST ary = E.parse();
    //if (ary instanceof ASTId) ary = Env.staticLookup((ASTId)ary);
    //// parse the probs, either a ASTSeries or an ASTSeq -> resulting in a Frame _ONLY_
    //AST seq = null;
    //// if is ASTSeries:
    //if (E.skipWS().peek() == '{') {
    //  String[] ps = E.xpeek('{').parseString('}').split(";");
    //  _probs = new double[ps.length];
    //  for (int i = 0; i < ps.length; ++i) {
    //    double v = Double.valueOf(ps[i]);
    //    if (v < 0 || v > 1) throw new  IllegalArgumentException("Quantile: probs must be in the range of [0, 1].");
    //    _probs[i] = v;
    //  }
    //
    //// else ASTSeq
    //} else seq = E.parse();
    //if (seq != null)
    //  if (seq instanceof ASTId) seq = Env.staticLookup((ASTId)seq);
    //// Get the na.rm
    //AST a = E._env.lookup((ASTId)E.skipWS().parse());
    //_narm = ((ASTNum)a).dbl() == 1;
    ////Get the type
    //try {
    //  _type = (int) ((ASTNum) E.skipWS().parse()).dbl();
    //} catch (ClassCastException e) {
    //  e.printStackTrace();
    //  throw new IllegalArgumentException("Argument `type` expected to be a number.");
    //}
    //// Finish the rest
    //ASTQtile res = (ASTQtile) clone();
    //res._asts = seq == null ? new AST[]{ary} : new AST[]{ary, seq}; // in reverse order so they appear correctly on the stack.
    //return res;
  }


  @Override protected void apply(Env env) {
    throw H2O.unimpl();
    //final Frame probs = _probs == null ? env.popAry() : null;
    //if (probs != null && probs.numCols() != 1) throw new IllegalArgumentException("Probs must be a single vector.");
    //
    //Frame x = env.popAry();
    //if (x.numCols() != 1) throw new IllegalArgumentException("Must specify a single column in quantile. Got: "+ x.numCols() + " columns.");
    //Vec xv  = x.anyVec();
    //if ( xv.isEnum() ) {
    //  throw new  IllegalArgumentException("Quantile: column type cannot be Categorical.");
    //}
    //
    //double p[];
    //
    //Vec pv = probs == null ? null : probs.anyVec();
    //if (pv != null) {
    //  p = new double[(int)pv.length()];
    //  for (int i = 0; i < pv.length(); i++) {
    //    if ((p[i] = pv.at((long) i)) < 0 || p[i] > 1)
    //      throw new IllegalArgumentException("Quantile: probs must be in the range of [0, 1].");
    //  }
    //} else p = _probs;
    //
    //String[] names = new String[p.length];
    //for (int i = 0; i < names.length; ++i) names[i] = Double.toString(p[i]) + "%";
    //
    //// create output vec
    //Vec res = Vec.makeZero(p.length);
    //
    //final int MAX_ITERATIONS = 16;
    //final int MAX_QBINS = 1000; // less uses less memory, can take more passes
    //final boolean MULTIPASS = true; // approx in 1 pass if false
    //// Type 7 matches R default
    //final int INTERPOLATION = _type; // 7 uses linear if quantile not exact on row. 2 uses mean.
    //
    //// a little obtuse because reusing first pass object, if p has multiple thresholds
    //// since it's always the same (always had same valStart/End seed = vec min/max
    //// some MULTIPASS conditionals needed if we were going to make this work for approx or exact
    //final Quantiles[] qbins1 = new Quantiles.BinningTask(MAX_QBINS, xv.min(), xv.max()).doAll(xv)._qbins;
    //for( int i=0; i<p.length; i++ ) {
    //  double quantile = p[i];
    //  // need to pass a different threshold now for each finishUp!
    //  qbins1[0].finishUp(xv, new double[]{quantile}, INTERPOLATION, MULTIPASS);
    //  if( qbins1[0]._done ) {
    //    res.set(i,qbins1[0]._pctile[0]);
    //  } else {
    //    // the 2-N map/reduces are here (with new start/ends. MULTIPASS is implied
    //    Quantiles[] qbinsM = new Quantiles.BinningTask(MAX_QBINS, qbins1[0]._newValStart, qbins1[0]._newValEnd).doAll(xv)._qbins;
    //    for( int iteration = 2; iteration <= MAX_ITERATIONS; iteration++ ) {
    //      qbinsM[0].finishUp(xv, new double[]{quantile}, INTERPOLATION, MULTIPASS);
    //      if( qbinsM[0]._done ) {
    //        res.set(i,qbinsM[0]._pctile[0]);
    //        break;
    //      }
    //      // the 2-N map/reduces are here (with new start/ends. MULTIPASS is implied
    //      qbinsM = new Quantiles.BinningTask(MAX_QBINS, qbinsM[0]._newValStart, qbinsM[0]._newValEnd).doAll(xv)._qbins;
    //    }
    //  }
    //}
    //res.chunkForChunkIdx(0).close(0,null);
    //res.postWrite(new Futures()).blockForPending();
    //Frame fr = new Frame(new String[]{"Quantiles"}, new Vec[]{res});
    //env.pushAry(fr);
  }
}


package hex.quantile;

import water.H2O;
import water.fvec.Frame;
import water.rapids.*;

class ASTMedian extends ASTReducerOp {
  ASTMedian() { super( 0 ); }
  @Override protected String opStr() { return "median"; }
  @Override protected ASTOp make() { return new ASTMedian(); }
  @Override protected ASTMedian parse_impl(Exec E) { return (ASTMedian)super.parse_impl(E); }
  @Override protected double op(double d0, double d1) { throw H2O.unimpl(); }
  @Override protected void apply(Env env) {
    Frame fr;
    try {
      fr = env.popAry();
    } catch (Exception e) {
      throw new IllegalArgumentException("`median` expects a single column from a Frame.");
    }
    if (fr.numCols() != 1)
      throw new IllegalArgumentException("`median` expects a single numeric column from a Frame.");

    if (!fr.anyVec().isNumeric())
      throw new IllegalArgumentException("`median` expects a single numeric column from a Frame.");

    throw H2O.unimpl();
    //Quantiles q = new Quantiles();
    //Quantiles[] qbins = new Quantiles.BinningTask(q._max_qbins, fr.anyVec().min(), fr.anyVec().max()).doAll(fr.anyVec())._qbins;
    //qbins[0].finishUp(fr.anyVec(), new double[]{0.5}, q._interpolation_type, true);
    //env.push(new ValNum(qbins[0]._pctile[0]));
  }
}

package hex.schemas;

import hex.quantile.QuantileModel;
import water.H2O;
import water.api.*;
import water.util.PojoUtils;

public class QuantileModelV2 extends ModelSchema<QuantileModel, QuantileModelV2, QuantileModel.QuantileParameters, QuantileV2.QuantileParametersV2, QuantileModel.QuantileOutput, QuantileModelV2.QuantileModelOutputV2> {

  public static final class QuantileModelOutputV2 extends ModelOutputSchema<QuantileModel.QuantileOutput, QuantileModelOutputV2> {
    @API(help = "Model parameters")
    QuantileV2.QuantileParametersV2 parameters;

    @API(help="Iterations executed")
    public double iterations;

    // Output fields
    @API(help="Quantiles")
    public double quantiles[/*N*/][/*Q*/]; // Our N columns, Q quantiles reported

    // TODO: gratuitous
    @Override public QuantileModel.QuantileOutput createImpl() {
      QuantileModel.QuantileOutput impl = new QuantileModel.QuantileOutput(null);
      PojoUtils.copyProperties(impl, this, PojoUtils.FieldNaming.DEST_HAS_UNDERSCORES);
      return impl;
    }

    // TODO: gratuitous
    @Override public QuantileModelOutputV2 fillFromImpl( QuantileModel.QuantileOutput impl) {
      PojoUtils.copyProperties(this, impl, PojoUtils.FieldNaming.ORIGIN_HAS_UNDERSCORES);
      return this;
    }
  } // QuantileModelOutputV2


  //==========================
  // Custom adapters go here

  // TOOD: I think we can implement the following two in ModelSchema, using reflection on the type parameters.
  public QuantileV2.QuantileParametersV2 createParametersSchema() { return new QuantileV2.QuantileParametersV2(); }
  public QuantileModelOutputV2 createOutputSchema() { return new QuantileModelOutputV2(); }

  // Version&Schema-specific filling into the impl
  @Override public QuantileModel createImpl() {
    throw H2O.unimpl();
  }

  // Version&Schema-specific filling from the impl
  @Override public QuantileModelV2 fillFromImpl( QuantileModel m ) {
    // TODO: if( !(h instanceof InspectHandler) ) throw H2O.unimpl();
    // TODO: InspectHandler ih = (InspectHandler)h;
    // TODO: QuantileModel kmm = ih._val.get();
    // TODO: iters = kmm._iters;
    return super.fillFromImpl(m);
  }
}

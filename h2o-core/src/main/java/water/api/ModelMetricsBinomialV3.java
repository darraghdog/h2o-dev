package water.api;

import hex.AUC2;
import hex.ModelMetricsBinomial;
import water.util.TwoDimTable;

public class ModelMetricsBinomialV3<I extends ModelMetricsBinomial, S extends ModelMetricsBinomialV3<I, S>> extends ModelMetricsBase<I,S> {
  @API(help="The standard deviation of the training response.", direction=API.Direction.OUTPUT)
    public double sigma; // Belongs in a mythical ModelMetricsSupervisedV3

  @API(help="The logarithmic loss for this scoring run.", direction=API.Direction.OUTPUT)
    public double logloss;

  @API(help="The AUC for this scoring run.", direction=API.Direction.OUTPUT)
    public double AUC;

  @API(help="The Gini score for this run.", direction=API.Direction.OUTPUT)
    public double Gini;

  @API(help="The Confusion Matrices for this scoring run; redundant but easier for the client to consume.", direction=API.Direction.OUTPUT)
  public long[][][] confusion_matrices;

  @API(help = "The Metrics for various thresholds.", direction = API.Direction.OUTPUT)
    public TwoDimTableBase thresholds_and_metric_scores;

  @API(help = "The Metrics for various criteria.", direction = API.Direction.OUTPUT)
    public TwoDimTableBase max_criteria_and_metric_scores;

  @Override
    public ModelMetricsBinomialV3 fillFromImpl(ModelMetricsBinomial modelMetrics) {
    super.fillFromImpl(modelMetrics);
    this.sigma = modelMetrics._sigma;
    this.logloss = modelMetrics._logloss;

    AUC2 auc = modelMetrics._auc;
    if (null != auc) {
      this.AUC  = auc._auc;
      this.Gini = auc._gini;
      
      // Fill TwoDimTable
      String[] thresholds = new String[auc._nBins];
      for( int i=0; i<auc._nBins; i++ )
        thresholds[i] = Double.toString(auc._ths[i]);
      AUC2.ThresholdCriterion crits[] = AUC2.ThresholdCriterion.VALUES;
      String[] colHeaders = new String[crits.length];
      String[] types      = new String[crits.length];
      String[] formats    = new String[crits.length];
      for( int i=0; i<crits.length; i++ ) {
        colHeaders[i] = crits[i].toString();
        types     [i] = crits[i]._isInt ? "long" : "double";
        formats   [i] = crits[i]._isInt ? "%d"   : "%f"    ;
      }
      TwoDimTable thresholdsByMetrics = new TwoDimTable("Thresholds x Metric Scores", null, thresholds, colHeaders, types, formats, "Thresholds" );
      confusion_matrices = new long[auc._nBins][2][2];

      for( int i=0; i<auc._nBins; i++ ) {
        for (int j = 0; j < crits.length; j++) {
          double d = crits[j].exec(auc, i); // Note: casts to Object are NOT redundant
          thresholdsByMetrics.set(i, j, crits[j]._isInt ? (Object) ((long) d) : d);
        }
        confusion_matrices[i] = auc.cm_for_threshold(i);
      }
      this.thresholds_and_metric_scores = new TwoDimTableV1().fillFromImpl(thresholdsByMetrics);
      
      // Fill TwoDimTable
      TwoDimTable maxMetrics = new TwoDimTable("Maximum Metric", null, colHeaders,
                                               new String[]{"Threshold","Metric","idx"},
                                               new String[]{"double",   "double","long"},
                                               new String[]{"%f",       "%f",    "%d"},
                                               "Metric" );
      for( int i=0; i<crits.length; i++ ) {
        int idx = crits[i].max_criterion_idx(auc);
        maxMetrics.set(i,0,idx==-1 ? Double.NaN : auc._ths[idx]);
        maxMetrics.set(i,1,idx==-1 ? Double.NaN : crits[i].exec(auc,idx));
        maxMetrics.set(i,2,idx);
      }
      
      this.max_criteria_and_metric_scores = new TwoDimTableV1().fillFromImpl(maxMetrics);
    }
    return this;
  }
}

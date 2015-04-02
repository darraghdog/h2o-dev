package water.rapids;

import water.Futures;
import water.Key;
import water.MRTask;
import water.fvec.*;

import java.util.ArrayList;

/**
* R's `apply`
*/
public class ASTApply extends ASTOp {
  protected static int _margin;  // 1 => work on rows; 2 => work on columns
  protected static String _fun;
  protected static AST[] _fun_args;
  static final String VARS[] = new String[]{ "", "ary", "MARGIN", "FUN", "..."};
  public ASTApply( ) { super(VARS); }
  @Override protected String opStr(){ return "apply";}
  @Override protected ASTOp make() {return new ASTApply();}
  @Override
  protected ASTApply parse_impl(Exec E) {
    AST ary = E.parse();
    if (ary instanceof ASTId) ary = Env.staticLookup((ASTId)ary);
    try {
      _margin = (int) ((ASTNum) E.skipWS().parse())._d;
    } catch (ClassCastException e) {
      throw new IllegalArgumentException("`MARGIN` must be either 1 or 2, it cannot be both.");
    }
    _fun = ((ASTId)E.skipWS().parse())._id;
    ArrayList<AST> fun_args = new ArrayList<>();
    while(E.skipWS().hasNext()) {
      fun_args.add(E.parse());
    }
    ASTApply res = (ASTApply)clone();
    res._asts = new AST[]{ary};
    if (fun_args.size() > 0) {
      _fun_args = fun_args.toArray(new AST[fun_args.size()]);
    } else {
      _fun_args = null;
    }
    return res;
  }
  @Override protected void apply(Env env) {
    String err="Result of function produced more than a single column!";
    // Peek everything from the stack
    final ASTOp op = ASTOp.get(_fun);
    Frame fr2 = null;  // results Frame
    Frame fr = env.popAry();
    env.addRef(fr);
    if( _margin == 2) {     // Work on columns?
      int ncols = fr.numCols();
      double[] row_result = new double[0];
      Vec[] vecs_result = new Vec[0];

      // Apply the function across columns
      // Types of results:
      //   A single row: Each col produces a single number result.
      //   A new array: Each column produces a new column
      //   If a new array, columns must be align'able.
      boolean isRow = false;

      // do the first column to determine the results type (isRow or not)
      Frame tmp = new Frame(new String[]{fr.names()[0]}, new Vec[]{fr.vecs()[0]});
      op.exec(env, new ASTFrame(tmp), _fun_args);
      if (env.isNum()) isRow = true;

      // if isRow, then append to row_result[]
      if (isRow) {
        row_result = new double[ncols];
        row_result[0] = env.popDbl();
      }

      // if !isRow, then append to vecs_result[]
      else {
        if (env.peekAry().numCols() != 1) throw new UnsupportedOperationException(err);
        vecs_result = new Vec[ncols];
        Frame v = env.popAry();
        vecs_result[0] = v.anyVec().makeCopy(null);
      }

      // loop over the columns and collect the results.
      // Appending to row_result or vecs_result accordingly
      for( int i=1; i<ncols; i++ ) {
        tmp = new Frame(new String[]{fr.names()[i]}, new Vec[]{fr.vecs()[i]});
        op.exec(env, new ASTFrame(tmp), _fun_args);
        if (isRow) row_result[i] = env.popDbl();
        else {
          if (env.peekAry().numCols() != 1) throw new UnsupportedOperationException(err);
          Frame v = env.popAry();
          vecs_result[i] = v.anyVec().makeCopy(null);
        }
      }

      // Create the results frame.
      if (isRow) {
        Futures fs = new Futures();
        Key key = Vec.VectorGroup.VG_LEN1.addVecs(1)[0];
        AppendableVec v = new AppendableVec(key);
        NewChunk chunk = new NewChunk(v, 0);
        for (double aRow_result : row_result) chunk.addNum(aRow_result);
        chunk.close(0, fs);
        Vec vec = v.close(fs);
        fs.blockForPending();
        fr2 = new Frame(vec);
      } else {
        fr2 = new Frame(fr.names(), vecs_result);
      }
    }
    if( _margin == 1) {      // Work on rows
      // apply on rows is essentially a map function
      // find out return type
      double[] rowin = new double[fr.vecs().length];
      for (int c = 0; c < rowin.length; c++) rowin[c] = fr.vecs()[c].at(0);
      final int outlen = op.map(env,rowin,null, _fun_args).length;
      final Env env0 = env;
      MRTask mrt = new MRTask() {
        @Override public void map(Chunk[] cs, NewChunk[] ncs) {
          double rowin [] = new double[cs.length];
          double rowout[] = new double[outlen];
          for (int row = 0; row < cs[0]._len; row++) {
            for (int c = 0; c < cs.length; c++) rowin[c] = cs[c].atd(row);
            rowout = op.map(env0, rowin, rowout, _fun_args);
            for (int c = 0; c < ncs.length; c++) ncs[c].addNum(rowout[c]);
          }
        }
      };
      String[] names = new String[outlen];
      for (int i = 0; i < names.length; i++) names[i] = "C"+(i+1);
      fr2 = mrt.doAll(outlen,fr).outputFrame(names, null);
    }
    else if (_margin != 1 && _margin != 2) throw new IllegalArgumentException("MARGIN limited to 1 (rows) or 2 (cols)");
    env.addRef(fr);
    env.pushAry(fr2);
  }
}

// --------------------------------------------------------------------------
// Same as "apply" but defaults to columns.
class ASTSApply extends ASTApply {
//  static final String VARS[] = new String[]{ "", "ary", "fcn", "..."};
  public ASTSApply( ) { super(); }
  @Override protected String opStr(){ return "sapply";}
  @Override protected ASTOp make() {return new ASTSApply();}
  @Override
  protected ASTSApply parse_impl(Exec E) {
    AST ary = E.parse();
    if (ary instanceof ASTId) ary = Env.staticLookup((ASTId)ary);
    _margin = 2;
    _fun = ((ASTId)E.skipWS().parse())._id;
    ArrayList<AST> fun_args = new ArrayList<>();
    while(E.skipWS().hasNext()) {
      fun_args.add(E.parse());
    }
    ASTSApply res = (ASTSApply)clone();
    res._asts = new AST[]{ary};
    if (fun_args.size() > 0) _fun_args = fun_args.toArray(new AST[fun_args.size()]);
    return res;
  }
  @Override protected void apply(Env env) {
    super.apply(env);
  }
}

// --------------------------------------------------------------------------
// unique(ary)
// Returns only the unique rows

//class ASTUnique extends ASTddply {
//  static final String VARS[] = new String[]{ "", "ary"};
//  ASTUnique( ) { super(VARS, new Type[]{Type.ARY, Type.ARY}); }
//  @Override protected String opStr(){ return "unique";}
//  @Override protected ASTOp make() {return new ASTUnique();}
//  @Override protected void apply(Env env, int argcnt, ASTApply apply) {
//    Thread cThr = Thread.currentThread();
//    Frame fr = env.peekAry();
//    int cols[] = new int[fr.numCols()];
//    for( int i=0; i<cols.length; i++ ) cols[i]=i;
//    ddplyPass1 p1 = new ddplyPass1( false, cols ).doAll(fr);
//    double dss[][] = new double[p1._groups.size()][];
//    int i=0;
//    for( Group g : p1._groups.keySet() )
//      dss[i++] = g._ds;
//    Frame res = FrameUtils.frame(fr._names,dss);
//    env.poppush(2,res,null);
//  }
//}

package water.rapids;

import water.DKV;
import water.Futures;
import water.H2O;
import water.Key;
import water.fvec.AppendableVec;
import water.fvec.Frame;
import water.fvec.NewChunk;
import water.fvec.Vec;

/**
 *  The ASTFunc Object
 *
 *  An ASTFunc pulls the function ast produced by the front-end and creates a reference to this function.
 *
 *  A function has a body (which may be empty), and a body is a list of statements.
 *
 *  Statements that are possible:
 *
 *  if statements
 *  else statements
 *  for statements
 *  while statements
 *  switch statement
 *  declarative statements
 *  operative statements
 *  return statements
 *
 *  The last statement of a function will return the result of that statement.
 *
 *  Some Rules:
 *  -----------
 *
 *  Every function defines a new Environment that inherits from the current one. Put another way, the calling scope
 *  provides the context for the function to be executed in. Environments can be captured with the `capture` call.
 *
 *  No function shall modify state in any of its parent environments (which includes the DKV store). A function may only
 *  return values to a parent scope.
 */


/**
 * Functions will always have all of their arguments fully specified (even if that means they are filled with null).
 * This means `nargs` arguments are always parsed.
 */
public class ASTFunc extends ASTFuncDef {
  private Env _e;  // a ref to the captures environment
  public ASTFunc() { super(); }
  public ASTFunc(String name, String[] arg_names, Env.SymbolTable table, ASTStatement body) {
    _name = name; _arg_names = arg_names; _table = table; _body = body;
  }
  AST[] _args;

  // (name args)
  @Override
  protected ASTFunc parse_impl(Exec E) {
    int nargs = _arg_names.length;
    AST[] args = new AST[nargs];
    for (int i = 0; i < nargs; ++i) {
      args[i] = E.skipWS().parse();
      if (args[i] instanceof ASTId) args[i] = Env.staticLookup((ASTId)args[i]);
    }

    ASTFunc res = (ASTFunc)clone();
    res._args = args;
    res._asts = _asts;
    return res;
  }

  @Override protected String opStr() { return _name; }

  @Override protected ASTOp make() { return new ASTFunc(); }

  @Override protected void apply(Env e) {
    Frame cleanme;
    Frame f;
    Env captured = e.capture();
    for (int i = 0; i < _args.length; ++i) {
      if (_args[i] instanceof ASTId) {
        ASTId a = (ASTId)_args[i];
        if (!a.isLookup()) throw new IllegalArgumentException("Function arguments must be lookups.");
        _args[i] = e.lookup(a);
      }
      if (!(_args[i] instanceof ASTNum) && !(_args[i] instanceof ASTString) && !(_args[i] instanceof ASTFrame) && !(_args[i] instanceof ASTNull)) {
        _args[i].treeWalk(e);
        _args[i] = e.pop2AST();
      }
      if (_args[i] instanceof ASTNum) _table.put(_arg_names[i], Env.NUM, _args[i].value());
      else if (_args[i] instanceof ASTString) _table.put(_arg_names[i], Env.STR, _args[i].value());
      else if (_args[i] instanceof ASTFrame) {
        // have a frame not in the DKV --> put in the DKV
        if (((ASTFrame)_args[i])._key == null) {
          cleanme = ((ASTFrame)_args[i])._fr;
          f = new Frame(Key.make(_arg_names[i]), cleanme.names(), cleanme.vecs());
          DKV.put(f._key, f); // block n put the key in the DKV
          _args[i] = new ASTFrame(f._key.toString());
          _table._frames.put(_arg_names[i], f);
        }
        _table.put(_arg_names[i], Env.LARY, _args[i].value());
        _table._frames.put(_arg_names[i], ((ASTFrame)_args[i])._fr);
      }
      else if (_args[i] instanceof ASTNull) _table.put(_arg_names[i], Env.STR, "null");
      else throw new IllegalArgumentException("Argument of type "+ _args[i].getClass()+" unsupported. Argument must be a String, number, Frame, or null.");
    }
    captured._local.copyOver(_table); // put the local table for the function into the _local table for the env
    _body.exec(captured);
//    captured.popScope();
    _e = captured;
  }

  // used by methods that pass their args to FUN (e.g. apply, sapply, ddply); i.e. args are not parsed here.
  @Override void exec(Env e, AST arg1, AST[] args) {
    _args = new AST[args == null ? 1 :1 + args.length];
    _args[0] = arg1;
    if (args != null) System.arraycopy(args, 0, _args, 1, args.length);
    apply(e);
  }

  double[] map(Env env, double[] in, double[] out, AST[] args) {
    Futures fs = new Futures();
    Vec[] vecs = new Vec[in.length];
    Key keys[] = Vec.VectorGroup.VG_LEN1.addVecs(vecs.length);
    for( int c = 0; c < vecs.length; c++ ) {
      AppendableVec vec = new AppendableVec(keys[c]);
      NewChunk chunk = new NewChunk(vec, 0);
      chunk.addNum(in[c]);
      chunk.close(0, fs);
      vecs[c] = vec.close(fs);
    }
    fs.blockForPending();
    Key local_key = Key.make();
    Frame fr = new Frame(local_key, null, vecs);
    _table._frames.put(local_key.toString(), fr); // push fr, since not in DKV, into the _local_frames -> must trash this frame at some point ... during popScope()

    // execute the function on the row
    exec(env, new ASTFrame(fr), args);

    // cleanup results and return
    if (env.isNum()) {
      if (out==null || out.length<1) out= new double[1];
      out[0] = env.popDbl();
    } else if (env.isAry()) {
      fr = env.popAry();
      if (fr.numCols() > 1 && fr.numRows() != 1) throw H2O.unimpl("Number of rows returned is > 1");
//      if (fr.numRows() > 1<<8) throw H2O.unimpl("Too many rows!");
      if (fr.numCols() > 1) {
        out = new double[fr.numCols()];
        for (int v = 0; v < fr.vecs().length; ++v) out[v] = fr.vecs()[v].at(0);
      } else {
        Vec vec = fr.anyVec();
        if (out == null || out.length < vec.length()) out = new double[(int) vec.length()];
        for (long i = 0; i < vec.length(); i++) out[(int) i] = vec.at(i);
      }
    } else {
      throw H2O.unimpl();
    }
//    env.cleanup(fr);
    return out;
  }

  @Override public StringBuilder toString( StringBuilder sb, int d ) {
    indent(sb,d).append(this).append(") {\n");
    _body.toString(sb,d+1).append("\n");
    return indent(sb,d).append("}");
  }

//  void trash() {
//    if (_e == null) return;
//    // wipe the _local_frames
//    Futures fs = new Futures();
//    for (Vec v: _e._refcnt.keySet()) {
//      fs = _e._parent.subVec(v, fs);
//    }
//    fs.blockForPending();
//  }
}

class ASTFuncDef extends ASTOp {
  protected static String _name;
  protected static String[] _arg_names;
  protected static Env.SymbolTable _table;
  protected ASTStatement _body;
  public ASTFuncDef() { super(null); }   // super(null) => _vars[] = null

  void parse_func(Exec E) {
    String name = E.parseID();
    _name = name;

    // parse the function args: these are just arg names -> will do _local.put(name, Env.NULL, null) (local ST put)
    Env.SymbolTable table = E._env.newTable(); // grab a new SymbolTable
    String[] args = E.skipWS().peek() == '{' ? E.xpeek('{').parseString('}').split(";") : null;
    for (int i = 0; i < args.length;++i) args[i] = args[i].replaceAll("\\s+","");
    _arg_names = args;
    if (args == null) table.put(null, Env.NULL, null);
    else for (String arg : args) table.put(arg, Env.NULL, null);
    _table = table;

    // parse the function body
    _body = new ASTStatement().parse_impl(E.skipWS());

    ASTFunc res = new ASTFunc(_name, _arg_names, _table, _body);
    res._asts = null;
    putUDF(res, name);  // not all nodes get this...
  }

  @Override protected String opStr() { return "def"; }

  @Override protected ASTOp make() { return new ASTFuncDef(); }

  @Override protected void apply(Env e) { throw H2O.fail(); }
}

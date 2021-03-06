package water.parser;

import java.io.IOException;
import java.io.InputStream;
import java.util.*;

import water.*;
import water.fvec.Vec;

/** A collection of utility classes for parsing.
 *
 *  Interfaces:
 *  DataIn - Manage bulk streaming input data to the parser.  Sometimes the data
 *           comes from parallel raw byte file reads, with speculative line
 *           starts.  Sometimes the data comes from an InputStream - probably a
 *           GZIP stream.
 *  DataOut- Interface for writing results of parsing, accumulating numbers and
 *           strings (enums) or handling invalid lines &amp; parse errors.
 *
 *  static classes:
 *  StreamData  - Class implementing DataIn from a Stream (probably a GZIP stream)
 *  InspectDataOut - Class implementing DataOut, on behalf of the GUI, for
 *                parsing &amp; previewing the first few lines &amp; columns of a file.
 */
public abstract class Parser extends Iced {
  static final byte CHAR_TAB = '\t';
  static final byte CHAR_CR = 13;
  static final byte CHAR_LF = 10;
  static final byte CHAR_SPACE = ' ';
  static final byte CHAR_DOUBLE_QUOTE = '"';
  static final byte CHAR_SINGLE_QUOTE = '\'';

  // State for the CSV & SVMLight Parser's FSAs
  protected static final byte SKIP_LINE = 0;
  protected static final byte EXPECT_COND_LF = 1;
  protected static final byte EOL = 2;
  protected static final byte TOKEN = 3;
  protected static final byte COND_QUOTED_TOKEN = 4;
  protected static final byte NUMBER = 5;
  protected static final byte NUMBER_SKIP = 6;
  protected static final byte NUMBER_SKIP_NO_DOT = 7;
  protected static final byte NUMBER_FRACTION = 8;
  protected static final byte NUMBER_EXP = 9;
  protected static final byte NUMBER_EXP_START = 11;
  protected static final byte NUMBER_END = 12;
  protected static final byte STRING = 13;
  protected static final byte COND_QUOTE = 14;
  protected static final byte SEPARATOR_OR_EOL = 15;
  protected static final byte WHITESPACE_BEFORE_TOKEN = 16;
  protected static final byte STRING_END = 17;
  protected static final byte COND_QUOTED_NUMBER_END = 18;
  protected static final byte POSSIBLE_EMPTY_LINE = 19;
  protected static final byte POSSIBLE_CURRENCY = 20;

  protected final byte CHAR_DECIMAL_SEP = '.';
  protected final byte CHAR_SEPARATOR;

  protected static final long LARGEST_DIGIT_NUMBER = Long.MAX_VALUE/10;
  protected static boolean isEOL(byte c) { return (c == CHAR_LF) || (c == CHAR_CR); }

  protected final ParseSetup _setup;
  Parser( ParseSetup setup ) { _setup = setup;  CHAR_SEPARATOR = setup._separator; }
  protected int fileHasHeader(byte[] bits, ParseSetup ps) { return ParseSetup.NO_HEADER; }

  // Parse this one Chunk (in parallel with other Chunks)
  abstract DataOut parseChunk(int cidx, final DataIn din, final DataOut dout);

  DataOut streamParse( final InputStream is, final DataOut dout) throws IOException {
    if( !_setup._parse_type._parallelParseSupported ) throw H2O.unimpl();
    StreamData din = new StreamData(is);
    int cidx=0;
    while( is.available() > 0 )
      parseChunk(cidx++, din, dout);
    parseChunk(cidx, din, dout);     // Parse the remaining partial 32K buffer
    return dout;
  }

  // ------------------------------------------------------------------------
  // Zipped file; no parallel decompression; decompress into local chunks,
  // parse local chunks; distribute chunks later.
  DataOut streamParseZip( final InputStream is, final StreamDataOut dout, InputStream bvs ) throws IOException {
    // All output into a fresh pile of NewChunks, one per column
    if( !_setup._parse_type._parallelParseSupported ) throw H2O.unimpl();
    StreamData din = new StreamData(is);
    int cidx=0;
    StreamDataOut nextChunk = dout;
    int zidx = bvs.read(null,0,0); // Back-channel read of chunk index
    assert zidx==1;
    while( is.available() > 0 ) {
      int xidx = bvs.read(null,0,0); // Back-channel read of chunk index
      if( xidx > zidx ) {  // Advanced chunk index of underlying ByteVec stream?
        zidx = xidx;       // Record advancing of chunk
        nextChunk.close(); // Match output chunks to input zipfile chunks
        if( dout != nextChunk ) dout.reduce(nextChunk);
        nextChunk = nextChunk.nextChunk();
      }
      parseChunk(cidx++, din, nextChunk);
    }
    parseChunk(cidx, din, nextChunk);     // Parse the remaining partial 32K buffer
    nextChunk.close();
    if( dout != nextChunk ) dout.reduce(nextChunk);
    return dout;
  }

  /** Manage bulk streaming input data to the parser.  Sometimes the data comes
   *  from parallel raw byte file reads, with speculative line starts.
   *  Sometimes the data comes from an InputStream - probably a GZIP stream.  */
  interface DataIn {
    // Get another chunk of byte data
    abstract byte[] getChunkData( int cidx );
    abstract int  getChunkDataStart( int cidx );
    abstract void setChunkDataStart( int cidx, int offset );
  }

  /** Interface for writing results of parsing, accumulating numbers and
   *  strings (enums) or handling invalid lines & parse errors.  */
  interface DataOut extends Freezable {
    void setColumnNames(String [] names);
    // Register a newLine from the parser
    void newLine();
    // True if already forced into a string column (skip number parsing)
    boolean isString(int colIdx);
    // Add a number column with given digits & exp
    void addNumCol(int colIdx, long number, int exp);
    // Add a number column with given digits & exp
    void addNumCol(int colIdx, double d);
    // An an invalid / missing entry
    void addInvalidCol(int colIdx);
    // Add a String column
    void addStrCol( int colIdx, ValueString str );
    // Final rolling back of partial line
    void rollbackLine();
    void invalidLine(String err);
  }

  interface StreamDataOut extends DataOut {
    StreamDataOut nextChunk();
    StreamDataOut reduce(StreamDataOut dout);
    StreamDataOut close();
    StreamDataOut close(Futures fs);
  }

  /** Class implementing DataIn from a Stream (probably a GZIP stream)
   *  Implements a classic double-buffer reader.
   */
  private static class StreamData implements Parser.DataIn {
    final transient InputStream _is;
    private byte[] _bits0 = new byte[64*1024];
    private byte[] _bits1 = new byte[64*1024];
    private int _cidx0=-1, _cidx1=-1; // Chunk #s
    private int _coff0=-1, _coff1=-1; // Last used byte in a chunk
    private StreamData(InputStream is){_is = is;}
    @Override public byte[] getChunkData(int cidx) {
      if( cidx == _cidx0 ) return _bits0;
      if( cidx == _cidx1 ) return _bits1;
      assert cidx==_cidx0+1 || cidx==_cidx1+1;
      byte[] bits = _cidx0<_cidx1 ? _bits0 : _bits1;
      if( _cidx0<_cidx1 ) { _cidx0 = cidx; _coff0 = -1; }
      else                { _cidx1 = cidx; _coff1 = -1; }
      // Read as much as the buffer will hold
      int off=0;
      try {
        while( off < bits.length ) {
          int len = _is.read(bits,off,bits.length-off);
          if( len == -1 ) break;
          off += len;
        }
        assert off == bits.length || _is.available() <= 0;
      } catch( IOException ioe ) {
        throw new RuntimeException(ioe);
      }
      if( off == bits.length ) return bits;
      // Final read is short; cache the short-read
      byte[] bits2 = (off == 0) ? null : Arrays.copyOf(bits,off);
      if( _cidx0==cidx ) _bits0 = bits2;
      else               _bits1 = bits2;
      return bits2;
    }
    @Override public int getChunkDataStart(int cidx) {
      if( _cidx0 == cidx ) return _coff0;
      if( _cidx1 == cidx ) return _coff1;
      return 0; 
    }
    @Override public void setChunkDataStart(int cidx, int offset) { 
      if( _cidx0 == cidx ) _coff0 = offset;
      if( _cidx1 == cidx ) _coff1 = offset;
    }
  }

  /** Class implementing DataOut, on behalf of the GUI, for parsing &amp;
   *  previewing the first few lines &amp; columns of a file.
   */
  protected static class InspectDataOut extends Iced implements DataOut {
    protected final static int MAX_PREVIEW_COLS  = 100;
    protected final static int MAX_PREVIEW_LINES = 10;
    protected int _nlines;
    protected int _ncols;
    protected int _invalidLines;
    private   String []   _colNames;
    protected String [][] _data = new String[MAX_PREVIEW_LINES][];
    protected boolean _dataIsColPreview = false;
    transient private HashSet<String> [] _domains;
    int [] _nnums;
    int [] _nstrings;
    int [] _ndates;
    int [] _nUUID;
    int [] _nzeros;
    int [] _nempty;
    transient ArrayList<String> _errors;

    protected InspectDataOut() {}
    protected InspectDataOut(int ncols) { setColumnCount(ncols); }

    String[] colNames() { return _colNames; }

    @Override public void setColumnNames(String[] names) {
      _colNames = names;
      _data[0] = names;
      ++_nlines;
      setColumnCount(names.length);
    }
    private void setColumnCount(int n) {
      // initialize
      if (_ncols == 0 && n > 0) {
        _ncols = n;
        _nzeros = new int[n];
        _nstrings = new int[n];
        _nUUID = new int[n];
        _ndates = new int[n];
        _nnums = new int[n];
        _nempty = new int[n];
        _domains = new HashSet[n];
        for(int i = 0; i < n; ++i)
          _domains[i] = new HashSet<String>();
        for(int i =0; i < MAX_PREVIEW_LINES; i++)
          _data[i] = new String[n];
      } /*else if (n > _ncols) { // resize
        _nzeros = Arrays.copyOf(_nzeros, n);
        _nstrings = Arrays.copyOf(_nstrings, n);
        _nUUID = Arrays.copyOf(_nUUID, n);
        _ndates = Arrays.copyOf(_ndates, n);
        _nnums = Arrays.copyOf(_nnums, n);
        _nempty = Arrays.copyOf(_nempty, n);
        _domains = Arrays.copyOf(_domains, n);
        for (int i=_ncols; i < n; i++)
          _domains[i] = new HashSet<String>();
        for(int i =0; i < MAX_PREVIEW_LINES; i++)
          _data[i] = Arrays.copyOf(_data[i], n);
        _ncols = n;
      }*/
    }
    @Override public void newLine() { ++_nlines; }
    @Override public boolean isString(int colIdx) { return false; }
    @Override public void addNumCol(int colIdx, long number, int exp) {
      if(colIdx < _ncols) {
        if (number == 0)
          ++_nzeros[colIdx];
        else
          ++_nnums[colIdx];
        if (_nlines < MAX_PREVIEW_LINES)
            _data[_nlines][colIdx] = Double.toString(number * water.util.PrettyPrint.pow10(exp));
      }
    }
    @Override public void addNumCol(int colIdx, double d) {
      if(colIdx < _ncols) {
        if (d == 0)
          ++_nzeros[colIdx];
        else
          ++_nnums[colIdx];
        if (_nlines < MAX_PREVIEW_LINES)
            _data[_nlines][colIdx] = Double.toString(d);
      }
    }
    @Override public void addInvalidCol(int colIdx) {
      if(colIdx < _ncols) {
        ++_nempty[colIdx];
        if (_nlines < MAX_PREVIEW_LINES)
            _data[_nlines][colIdx] = "NA";
      }
    }
    @Override public void addStrCol(int colIdx, ValueString str) {
      if(colIdx < _ncols) {
        // Check for time
        if (ParseTime.isDateTime(str)) {
          ++_ndates[colIdx];
          return;
        }

        //Check for UUID
        if(ParseTime.isUUID(str)) {
          ++_nUUID[colIdx];
          return;
        }

        //Add string to domains list for later determining string, NA, or enum
        ++_nstrings[colIdx];
        _domains[colIdx].add(str.toString());

        if (_nlines < MAX_PREVIEW_LINES)
            _data[_nlines][colIdx] = str.toString();
      }
    }

    @Override public void rollbackLine() {--_nlines;}
    @Override public void invalidLine(String err) {
      ++_invalidLines;
      if( _errors == null ) _errors = new ArrayList<>();
      if( _errors.size() < 10 )
        _errors.add("Error at line: " + _nlines + ", reason: " + err);
    }
    String[] errors() { return _errors == null ? null : _errors.toArray(new String[_errors.size()]); }

    public byte[] guessTypes() {
      byte[] types = new byte[_ncols];
      for (int i = 0; i < _ncols; ++i) {
        int nonemptyLines = _nlines - _nempty[i] - 1; //During guess, some columns may be shorted one line based on 4M boundary

        //Very redundant tests, but clearer and not speed critical
        
        // is it clearly numeric?
        if ((_nnums[i] + _nzeros[i]) >= _ndates[i]
                && (_nnums[i] + _nzeros[i]) >= _nUUID[i]
                && _nnums[i] >= _nstrings[i]) { // 0s can be an NA among enums, ignore
          types[i] = Vec.T_NUM;
          continue;
        }

        // with NA, but likely numeric
        if (_domains[i].size() <= 1
                && (_nnums[i] + _nstrings[i] + _nzeros[i]) > _ndates[i] + _nUUID[i]) {
          types[i] = Vec.T_NUM;
          continue;
        }

        // Datetime
        if (_ndates[i] > _nUUID[i]
                && _ndates[i] > (_nnums[i] + _nzeros[i])
                && (_ndates[i] > _nstrings[i] || _domains[i].size() <= 1)) {
          types[i] = Vec.T_TIME;
          continue;
        }

        // UUID
        if (_nUUID[i] > _ndates[i]
                && _nUUID[i] > (_nnums[i] + _nzeros[i])
                && (_nUUID[i] > _nstrings[i] || _domains[i].size() <= 1)) {
          types[i] = Vec.T_UUID;
          continue;
        }

        // Strings, almost no dups
        if (_nstrings[i] > _ndates[i]
                && _nstrings[i] > _nUUID[i]
                && _nstrings[i] > (_nnums[i] + _nzeros[i])
                && _domains[i].size() >= 0.95 * _nstrings[i]) {
          types[i] = Vec.T_STR;
          continue;
        }

        // Enum or string?
        // Enum with 0s for NAs
        if(_nzeros[i] > 0
                && ((_nzeros[i] + _nstrings[i]) >= nonemptyLines) //just strings and zeros for NA (thus no empty lines)
                && (_domains[i].size() <= 0.95 * _nstrings[i]) ) { // not all unique strings
          types[i] = Vec.T_ENUM;
          continue;
        }
        // Enum mixed with numbers
        if(_nstrings[i] >= (_nnums[i]+_nzeros[i]) // mostly strings
                && (_domains[i].size() <= 0.95 * _nstrings[i]) ) { // but not all unique
          types[i] = Vec.T_ENUM;
          continue;
        }

        // All guesses failed
        types[i] = Vec.T_BAD;
      }
      return types;
    }

    public String[] guessNAStrings(byte[] types) {
      //For now just catch 0's as NA in Enums
      String[] na_strings = new String[_ncols];
      for (int i = 0; i < _ncols; ++i) {
        int nonemptyLines = _nlines - _nempty[i] - 1; //During guess, some columns may be shorted one line (based on 4M boundary)
        if (types[i] == Vec.T_ENUM
                && _nzeros[i] > 0
                && ((_nzeros[i] + _nstrings[i]) >= nonemptyLines) //just strings and zeros for NA (thus no empty lines)
                && (_domains[i].size() <= 0.95 * _nstrings[i])) { // not all unique strings
          na_strings[i] = "0";
        }
      }
      return na_strings;
    }
  }
}

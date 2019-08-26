//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *                  Marco Esters - Duke University 2019                    *
// *                                                                         *
//****************************************************************************

#include "aflowlib.h"

#define _DB_WEB_DEBUG_ false

using std::string;
using std::vector;

static const string _WEBENTRY_TABLE_ = "WEBENTRY";

namespace aflowlib {

//getAflowlibEntryWeb/////////////////////////////////////////////////////////
// Wrapper function to return the webentry in JSON
void getAflowlibEntryWeb(const string& db_file, const string& id, std::ostream& oss) {
  AflowDB db(db_file);
  oss << getEntryJSON(db, id, true);
}

//getEntryJSON////////////////////////////////////////////////////////////////
// Formats the data of an entry as a JSON object. The Boolean html determines
// whether the entries should be formatted using HTML tags or if they should
// be returned in their raw form.
string getEntryJSON(AflowDB& db, string id, bool html) {
  // Determine the format of the input ID.
  string search_by;
  if (aurostd::substring2bool(id, "aflow:")) {
    search_by = "auid";
    id = aurostd::tolower(id);
  } else if (aurostd::substring2bool(id, "aflowlib")) {
    search_by = "aurl";
  } else if (aurostd::substring2bool(id, "icsd:")) {
    search_by = "icsd";
    id = aurostd::RemoveSubStringFirst(id, "icsd:");
  } else {  // If neither, assume it is a title
    search_by = "title";
  }
  
  db.transaction(true);

  // Prepare temporary table containing only the desired entry.
  // This makes retrieving values much quicker.
  vector<string> tables = db.getTableSubset("schema");
  uint ntables = tables.size();

  string as, as_where;
  if (search_by == "icsd") as_where = "prototype LIKE '%_" + id + "'";
  else as_where = search_by + "='" + id + "'";

  if (search_by == "auid") {
    as = "SELECT * FROM auid_" + id.substr(6, 2) + " WHERE " + as_where;
  } else {
    for (uint t = 0; t < ntables; t++) {
      as += "SELECT * FROM " + tables[t] + " WHERE " + as_where;
      if (t < ntables - 1) as += " UNION ALL ";
    }
  }
  as += "LIMIT 1";
  db.createTempTableAs(_WEBENTRY_TABLE_, as);

  // Return empty string if no entry found
  if (aurostd::string2utype<int>(db.getProperty("COUNT", _WEBENTRY_TABLE_, "*")) == 0) {
    db.dropTempTable(_WEBENTRY_TABLE_);
    db.transaction(false);
    return "";
  }

  // Get all properties with production status. Apart from the property name,
  // the data type (for proper JSON formatting) and the formatter function are
  // also required.
  string where = "status='production'";
  vector<string> cols;
  aurostd::string2tokens("name,type,function", cols, ",");
  vector<vector<string> > schema = db.getSetMultiCol(_DB_SCHEMA_TABLE_, cols, false, where);
  uint nschema = schema.size();

  // The files list is important for links and images
  string files = db.getValue(_WEBENTRY_TABLE_, "files");

  std::stringstream json;
  string tab = "  ";
  string value;

  // Finally build the JSON output
  json << "{" << std::endl;
  for (uint s = 0; s < nschema; s++) {
    json << tab << "\"" << schema[s][0] << "\": ";
    // Links and images are handled differently than data
    if ((schema[s][2] == "link") || (schema[s][2] == "image")) {
      json << formatHREFToJSON(schema[s][0], files, db);
    } else {
      if (schema[s][0] != "icsd_number") {  // ICSD number is not stored in the database
        value = db.getValue(_WEBENTRY_TABLE_, schema[s][0]);
      }
      if (html || (schema[s][0] == "icsd_number")) {
        value = formatValueToHTML(value, schema[s][2], db);
      }
      if (value.empty()) {
        json << "null";
      } else if (schema[s][1] == "string") {
        json << "\"" << value << "\"";
      } else {
        json << value;
      }
    }
    if (s < nschema - 1) json << ",";
    json << std::endl;
  }
  json << "}" << std::endl;

  // Replace some placeholders from the schema
  string prototype = db.getValue(_WEBENTRY_TABLE_, "prototype");
  aurostd::StringStreamSubst(json, "$prototype", prototype);
  string aurl = db.getValue(_WEBENTRY_TABLE_, "aurl");
  aurostd::StringSubst(aurl, ":AFLOWDATA", "/AFLOWDATA");  // aurl has : instead of / after duke.edu
  aurostd::StringStreamSubst(json, "$aurl", aurl);
  // lattice_orig has to be substituted first!
  string lattice_orig = db.getValue(_WEBENTRY_TABLE_, "Bravais_lattice_orig");
  aurostd::StringStreamSubst(json, "$lattice_orig", lattice_orig);
  string lattice = db.getValue(_WEBENTRY_TABLE_, "Bravais_lattice_relax");
  aurostd::StringStreamSubst(json, "$lattice", lattice);

  // Format some special cases to HTML
  if (html) {
    aurostd::StringStreamSubst(json, "\\\\Gamma", "&Gamma;");
    aurostd::StringStreamSubst(json, "2_x_", "2 &times; ");
  }

  db.dropTempTable(_WEBENTRY_TABLE_);
  db.transaction(false);
  return json.str();
}

}  // namespace aflowlib

/******************************** FORMATTERS ********************************/

namespace aflowlib {

//formatHREFToJSON////////////////////////////////////////////////////////////
// Returns a JSON object that contains the href/src and the link text/alt text
// of a link/image.
string formatHREFToJSON(const string& name, const string& files, AflowDB& db) {
  string href, linktext, filename, where;

  where = "name='" + name + "'";
  href = db.getValue(_DB_SCHEMA_TABLE_, "href", where);
  linktext = "\"" + db.getValue(_DB_SCHEMA_TABLE_, "linktext", where) + "\"";

  // Check if the linked file exists
  filename = aurostd::RemoveSubStringFirst(href, "$aurl/");
  // Replace $prototype placeholder (e.g. for Bader)
  if (aurostd::substring2bool(filename, "$prototype")) {
    string prototype = db.getValue(_WEBENTRY_TABLE_, "prototype");
    aurostd::StringSubst(filename, "$prototype", prototype);
  }
  // bz_image is not inside the files list, but the image is always on the server
  if (aurostd::substring2bool(name, "bz_image") || aurostd::substring2bool(files, filename)) {
    href = "\"" + href + "\"";
  } else {
    href = "null";
  }

  return "{\"href\": " + href + ", \"linktext\": " + linktext + "}";
}

//formatValueToHTML///////////////////////////////////////////////////////////
// formats 
string formatValueToHTML(const string& value, const string& function, AflowDB& db) {
  if (function.empty()) return value;  // No formatting function, so return as is
  if (function == "compound") return formatCompound(value);
  if (function == "crystallography") return formatCrystallography(value);
  if (function == "icsd") return formatICSD(db);
  if (function == "lattice") return formatLattice(value, false); 
  if (function == "lattice_reciprocal") return formatLattice(value, true);
  if (function == "ldau") return formatLDAU(value);
  if (function == "spacegroup") return formatSpaceGroup(value);
  if (function == "species_resolved") return formatSpeciesResolvedProperties(value, db);
  if (function == "atom_resolved") return formatAtomResolvedProperties(value, db);
  // If the code arrives here, the function has not been implemented.
  // Return as is, but display a warning.
  std::cerr << "WARNING: function " << function << " not implemented." << std::endl;
  return value;
}

//formatCompound//////////////////////////////////////////////////////////////
// Takes a chemical formula and returns it with the appropriate subscripts.
string formatCompound(const string& label) {
  string html;
  bool sub = false;  // <sub> open?
  for (uint i = 0; i < label.size(); i++) {
    if (isdigit(label[i]) && !sub) {
      html += "<sub>";
      sub = true;
    } else if (isalpha(label[i]) && sub) {
      html += "</sub>";
      sub = false;
    }
    html += label[i];
  }
  // Close tags that may still be open
  if (sub) html += "</sub>";
  return html;
}

//formatCrystallography///////////////////////////////////////////////////////
// Applies formatting according to conventions in crystallography for
// space groups and point groups (Hermann-Mauguin and Schoenflies).
string formatCrystallography(const string& label) {
  string html;
  bool it = false;   // <i> open?
  bool sub = false;  // <sub> open?
  bool ol = false;   // overline open?
  for (uint i = 0; i < label.size(); i++) {
    if (isalpha(label[i])) {
      if (!it) {
        html += "<i>";
        it = true;
      }
      html += label[i];
    } else {
      if (it) {
        html += "</i>";
        it = false;
      }
      // This assumes that everything after _ is subscript unless it is
      // enclosed by {}. This is the current convention in AFLOW.
      if (label[i] == '_') {
        html += "<sub>";
        sub = true;
        if (label[i+1] == '{') i++;  // Skip { character
      } else if (label[i] == '}') {
        html += "</sub>";
        sub = false;
      } else if (label[i] == '-') {
        html += "<span style=\\\"text-decoration:overline;\\\">";
        ol = true;
      } else {
        html += label[i];
        if (isdigit(label[i]) && ol) {  // Close overline
          html += "</span>";
          ol = false;
        }
      }
    }
  }
  // Close tags that may still be open.
  if (it) html += "</i>";
  if (sub) html += "</sub>";
  return html;
}

//formatICSD//////////////////////////////////////////////////////////////////
// Returns the ICSD number using the prototype.
string formatICSD(AflowDB& db) {
  string catalog = db.getValue(_WEBENTRY_TABLE_, "catalog");
  if (catalog != "ICSD") return "";
  string prototype = db.getValue(_WEBENTRY_TABLE_, "prototype");
  vector<string> tokens;
  aurostd::string2tokens(prototype, tokens, "_");
  return tokens.back();
}

//formatLattice///////////////////////////////////////////////////////////////
// Takes an array of (reciprocal) lattice parameters and returns a formatted
// HTML string.
string formatLattice(const string& lattice_string, bool reciprocal) {
  if (lattice_string.empty()) return "";
  // Convert the geometry array into a vector
  vector<string> lattice = getVectorFromArrayString(lattice_string);

  // Prepare the labels
  vector<string> axes, angles;
  if (reciprocal) {
    aurostd::string2tokens("a<sup>*</sup>,b<sup>*</sup>,c<sup>*</sup>", axes, ",");
  } else {
    aurostd::string2tokens("a,b,c", axes, ",");
  }
  aurostd::string2tokens("&alpha;,&beta;,&gamma;", angles, ",");

  // Write lattice information
  string html = "\"";
  for (int i = 0; i < 3; i++) {
    html += "<i>" + axes[i] + "</i> = " + lattice[i] + " &Aring;";
    if (reciprocal) html += "<sup>-1</sup>";
    if (i < 2) html += "; ";
  }
  if (!reciprocal) {
    double ca_ratio = aurostd::string2utype<double>(lattice[2])/aurostd::string2utype<double>(lattice[0]);
    html += "; <i>c</i>/<i>a</i> = " + aurostd::utype2string<double>(ca_ratio) + ";";
  }
  for (int i = 0; i < 3; i++) {
    html += "<i>" + angles[i] + "</i> = " + lattice[i+3] + "&deg;";
    if (i < 2) html += "; ";
  }
  html += "\"";
  return html;
}

//formatLDAU//////////////////////////////////////////////////////////////////
// Returns the LDAU method from the type number (VASP convention).
string formatLDAU(const string& ldau_type) {
  int ldau;
  if (ldau_type.empty()) ldau = 0;
  else ldau = aurostd::string2utype<int>(ldau_type);
  if (ldau == 0) return "\"no LDAU\"";
  if (ldau == 1) return "\"1 (Liechtenstein LSDAU)\"";
  if (ldau == 2) return "\"2 (Dudarev)\"";
  if (ldau == 4) return "\"4 (Liechtenstein LDAU)\"";
  return "";
}

//formatSpaceGroup////////////////////////////////////////////////////////////
// Returns the spacegroup number along with its symbol.
string formatSpaceGroup(string number) {
  if (number.empty()) return "";
  number = aurostd::utype2string<int>(aurostd::string2utype<int>(number));
  string name = GetSpaceGroupName(aurostd::string2utype<int>(number));
  return string("\"" + formatCrystallography(name) + " (# " + number + ")\"");
}

//formatSpeciesResolvedProperties/////////////////////////////////////////////
// Wrapper function for AtomResolvedProperties where each specie has only
// one property.
string formatSpeciesResolvedProperties(const string& properties, AflowDB& db) {
  if (properties.empty()) return "";
  string val;
  val = db.getValue(_WEBENTRY_TABLE_, "species");
  vector<string> species = getVectorFromArrayString(val);
  vector<int> composition(species.size(), 1);
  return formatAtomResolvedProperties(species, composition, properties, false);
}

//formatAtomResolvedProperties////////////////////////////////////////////////
// This function takes a flat array and assigns the properties in that array
// to the atoms in the structure.

// This wrapper function prepares the composition vector
string formatAtomResolvedProperties(string properties, AflowDB& db) {
  if (properties.empty()) return "";
  string val;
  val = db.getValue(_WEBENTRY_TABLE_, "species");
  vector<string> species = getVectorFromArrayString(val);
  vector<int> composition(species.size());
  // If the properties are a 2D-array, the composition is implicitly given
  // by the array. So, extract the array dimensions and flatten it.
  if (aurostd::substring2bool(properties, "[[")) {
    vector<string> tokens, tokens2;
    properties = aurostd::RemoveSubStringFirst(properties, "[[");
    properties = aurostd::RemoveSubString(properties, "]");
    aurostd::string2tokens(properties, tokens, "[");
    for (uint i = 0; i < tokens.size(); i++) {
      aurostd::string2tokens(tokens[i], tokens2, ",");
      composition[i] = (int) tokens2.size();  // Trailing commas are eliminated by string2tokens
    }
    properties = "[" + aurostd::RemoveSubString(properties, "[") + "]";
  // Else get the composition from the database
  } else {
    val = db.getValue(_WEBENTRY_TABLE_, "composition");
    composition = aurostd::vectorstring2vectorint(getVectorFromArrayString(val));
  }
  return formatAtomResolvedProperties(species, composition, properties, true);
}

// This function returns an array of JSON objects that store the species and
// the property as an array or a single value.
string formatAtomResolvedProperties(const vector<string>& species,
                                    const vector<int>& composition, 
                                    const string& properties,
                                    bool array) {
  if (properties.empty()) return "";
  std::stringstream out;
  int iatom = 0;
  vector<string> prop = getVectorFromArrayString(properties);
  out << "[";
  for (uint s = 0; s < species.size(); s++) {
    out << "{\"species\": " << species[s] << ", \"value\": ";
    if (array) out << "[";
    for (int i = 0; i < composition[s]; i++) {
      out << prop[iatom];
      if (i < composition[s] - 1) out << ", ";
      iatom++;
    }
    if (array) out << "]";
    out << "}";
    if (s < species.size() - 1) out << ", ";
  }
  out << "]";
  return out.str();
}

//getVectorFromArrayString////////////////////////////////////////////////////
// Takes string representing a 1D-array and returns it as a vector of strings.
vector<string> getVectorFromArrayString(string arr) {
  vector<string> varray;
  arr = aurostd::RemoveSubString(arr, "[");
  arr = aurostd::RemoveSubString(arr, "]");
  aurostd::string2tokens(arr, varray, ",");
  return varray;
}

}  // namespace aflowlib

//****************************************************************************
// *                                                                         *
// *           Aflow STEFANO CURTAROLO - Duke University 2003-2019           *
// *                  Marco Esters - Duke University 2019                    *
// *                                                                         *
//****************************************************************************

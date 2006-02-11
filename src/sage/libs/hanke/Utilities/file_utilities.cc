

/////////////////////////////
// Checks if a file exists //
/////////////////////////////

bool FileExists(char* filename) {

  ifstream file;
  file.open(filename);

  // Check if the file exists
  if (file.is_open() == true) {
    file.close();
    return true;
  }

  return false;

}


bool FileExists(const char* & filename) {

  ifstream file;
  file.open(filename);

  // Check if the file exists
  if (file.is_open() == true) {
    file.close();
    return true;
  }

  return false;

}


bool FileExists(string filename) {

  ifstream file;
  file.open(filename.c_str());

  // Check if the file exists
  if (file.is_open() == true) {
    file.close();
    return true;
  }

  return false;

}




//////////////////////////////////
// Checks if a directory exists //   <== NOTE: This is a copy of FileExists()!
//////////////////////////////////

bool DirectoryExists(char* filename) {

  ifstream file;
  file.open(filename);

  // Check if the file exists
  if (file.is_open() == true) {
    file.close();
    return true;
  }

  return false;

}


bool DirectoryExists(const char* & filename) {

  ifstream file;
  file.open(filename);

  // Check if the file exists
  if (file.is_open() == true) {
    file.close();
    return true;
  }

  return false;

}


bool DirectoryExists(string filename) {

  ifstream file;
  file.open(filename.c_str());

  // Check if the file exists
  if (file.is_open() == true) {
    file.close();
    return true;
  }

  return false;

}




///////////////////////////////////////////////////
// Returns the absolute path to a file/directory //
///////////////////////////////////////////////////

string GetAbsolutePath(string filename) {

  // Check if the given path is already absolute
  if (filename[0] == '/') {
    //    string newfilename(filename);
    //    return newfilename;

    return string(filename);
  }

  // If not, then either find the path to "~/" or to "./"
  string newfilename(filename);
  string filepath;
  if ((filename[0] == '~') && (filename[1] == '/')) {
    filepath = string(getenv("HOME"));           // Find the home directory path
    newfilename.erase(0,2);                      // Forget the first two characters
  }
  else if ((filename[0] == '.') && (filename[1] == '/')) {
    filepath = string(getenv("PWD"));            // Find the present working directory path
  }


  // Concatenate this with the given filename
  string absolutepath;
  absolutepath = filepath + "/" + newfilename;

  // Return the absolute path
  return absolutepath;

}




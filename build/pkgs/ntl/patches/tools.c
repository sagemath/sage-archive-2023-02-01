
#include <NTL/tools.h>

#include <ctype.h>

#include <NTL/new.h>

NTL_START_IMPL


/*
   The following code differs from vanilla NTL 5.4.2.

   We add a SetErrorCallbackFunction(). This sets a global callback function _function_,
   which gets called with parameter _context_ and an error message string whenever Error()
   gets called.

   Note that if the custom error handler *returns*, then NTL will dump the error message
   back to stderr and abort() as it habitually does.

   -- David Harvey (2008-04-12)
*/

void (*ErrorCallbackFunction)(const char*, void*) = NULL;
void *ErrorCallbackContext = NULL;


void SetErrorCallbackFunction(void (*function)(const char*, void*), void *context)
{
   ErrorCallbackFunction = function;
   ErrorCallbackContext = context;
}


void Error(const char *s)
{
   if (ErrorCallbackFunction != NULL)
      ErrorCallbackFunction(s, ErrorCallbackContext);

   cerr << s << "\n";
   abort();
}


// The following implementation of CharToIntVal is completely portable.

long CharToIntVal(long a)
{
   switch (a) {
      case '0': return 0;
      case '1': return 1;
      case '2': return 2;
      case '3': return 3;
      case '4': return 4;
      case '5': return 5;
      case '6': return 6;
      case '7': return 7;
      case '8': return 8;
      case '9': return 9;

      case 'A': return 10;
      case 'B': return 11;
      case 'C': return 12;
      case 'D': return 13;
      case 'E': return 14;
      case 'F': return 15;

      case 'a': return 10;
      case 'b': return 11;
      case 'c': return 12;
      case 'd': return 13;
      case 'e': return 14;
      case 'f': return 15;

      default:  return -1;
   }
}

// The following implementation of IntValToChar is completely portable.

char IntValToChar(long a)
{
   switch (a) {
      case 0: return '0';
      case 1: return '1';
      case 2: return '2';
      case 3: return '3';
      case 4: return '4';
      case 5: return '5';
      case 6: return '6';
      case 7: return '7';
      case 8: return '8';
      case 9: return '9';

      case 10: return 'a';
      case 11: return 'b';
      case 12: return 'c';
      case 13: return 'd';
      case 14: return 'e';
      case 15: return 'f';

      default: Error("IntValToChar: bad arg");
   }

   return 0;  // to supress warnings
}


long IsWhiteSpace(long a)
{
   if (a > NTL_MAX_INT || a < NTL_MIN_INT)
      return 0;

   int b = (int) a;

   if (isspace(b))
      return 1;
   else
      return 0;
}

long SkipWhiteSpace(istream& s)
{
   long c;

   c = s.peek();
   while (IsWhiteSpace(c)) {
      s.get();
      c = s.peek();
   }

   if (c == EOF)
      return 0;
   else
      return 1;
}


void PrintTime(ostream& s, double t)
{
   long hh, mm, ss;

   ss = long(t + 0.5);

   hh = ss/3600;
   ss = ss - hh*3600;
   mm = ss/60;
   ss = ss - mm*60;

   if (hh > 0)
      s << hh << ":";

   if (hh > 0 || mm > 0) {
      if (hh > 0 && mm < 10) s << "0";
      s << mm << ":";
   }

   if ((hh > 0 || mm > 0) && ss < 10) s << "0";
   s << ss;
}

NTL_END_IMPL

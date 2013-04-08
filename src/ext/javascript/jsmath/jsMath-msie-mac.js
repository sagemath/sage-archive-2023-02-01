/*
 *  jsMath-msie-mac.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file makes changes needed for use with MSIE on the Mac.
 *
 *  ---------------------------------------------------------------------
 *
 *  Copyright 2004-2006 by Davide P. Cervone
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */



/*
 *  MSIE crashes if it changes the page too quickly, so we add a
 *  delay between processing math entries.  Unfortunately, this really
 *  slows down math in MSIE on the mac.
 */

jsMath.Add(jsMath,{

  msieProcess: jsMath.Process,
  msieProcessBeforeShowing: jsMath.ProcessBeforeShowing,

  Process: function () {
    // we need to delay a bit before starting to process the page
    //   in order to avoid an MSIE display bug
    jsMath.Message.Set("Processing Math: 0%");
    setTimeout('jsMath.msieProcess()',jsMath.Browser.delay);
  },

  ProcessBeforeShowing: function () {
    // we need to delay a bit before starting to process the page
    //   in order to avoid an MSIE display bug
    setTimeout('jsMath.msieProcessBeforeShowing()',5*jsMath.Browser.delay);
  }

});

jsMath.Browser.delay = 75;  // hope this is enough of a delay!

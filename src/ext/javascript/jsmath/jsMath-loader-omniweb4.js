/*
 *  jsMath-loader-omniweb4.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file synchronizes the jsMath-loader.html file with
 *  the actual loading of the source javascript file.
 *  OmniWeb 4 has a serious bug where the loader file is run
 *  several times (and out of sequence), which plays havoc
 *  with the Start() and End() calls.
 *
 *  ---------------------------------------------------------------------
 *
 *  Copyright 2006 by Davide P. Cervone
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

if (window.jsMathAutoload) {
  jsMath.Autoload.Script.endLoad();
} else {
  if (!window.phase2) {
    jsMath.Script.Start();
    window.phase2 = 1;
  } else {
    jsMath.Script.End();
    jsMath.Script.endLoad();
  }
}
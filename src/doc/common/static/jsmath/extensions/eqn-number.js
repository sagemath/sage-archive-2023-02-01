/*
 *  extensions/eqn-number.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file causes jsMath to add equation numbers to displayed
 *  equations.  These are displayed at the right, but the styles can
 *  be controlled through the jsMath.EqnNumber object.  Equations
 *  are numbered if they include a \label{xxx} call, and the macro
 *  \ref{xxx} can be used to refer to the equation number elsewhere
 *  in the document (it must appear by itself in a math formula,
 *  e.g., $\ref{xxx}$).  The "label-ref" CSS style can be used to
 *  style the references.
 *
 *  If jsMath.EqnNumber.autonumber is set to 1, then ALL displayed
 *  equations will be numberd.  Use the \nolabel macro to prevent
 *  equation numbering on an equation.
 *
 *  You can activate eqn-numbering by calling
 *
 *    jsMath.Extension.Require('eqn-number');
 *
 *  once jsMath.js has been loaded, or by adding "extensions/eqn-number.js"
 *  to the loadFiles array in jsMath/easy/load.js.
 *
 *  ---------------------------------------------------------------------
 *
 *  Copyright 2008 by Davide P. Cervone
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

/********************************************************************/

if (jsMath.EqnNumber) {jsMath.EqnNumber_old = jsMath.EqnNumber}

jsMath.EqnNumber = {

  styles: {
    '.jsMath_displayBox, .tex2math_div': {position: 'relative'},
    '.jsMath_number': {
        position: 'absolute',
        right: '2em', top: '50%', 'margin-top': '-.5em',
        height: 'auto', width: 'auto'
    },
    '.jsMath_ref': {'text-decoration': 'none'}
  },

  autonumber: 0,  // set to 1 to have ALL equations numbered

  number: 0,
  format: function (n) {return n},
  formatLabel: function (n) {return '<A NAME="eqn-'+n+'">('+n+')</A>'},
  formatRef: function (n) {return '(<A CLASS="jsMath_ref" HREF="#eqn-'+n+'">'+n+'</A>)'},

  _label: null,   // flag set when \label{x} is used
  _labels: {},    // stores label-name => label-value pairs
  _refs: {},      // stores elements referring to undefined labels
  _nolabel: 0,    // set by \nolabel

  nextNumber: function () {
    var ref = this.format(++this.number);
    if (this._label) {
      this._labels[this._label] = ref;
      if (this._refs[this._label]) this.fixRefs(this._label);
    }
    return this.formatLabel(ref);
  },

  isRef: function (element) {
    var tex = element.innerHTML;
    var result = tex.match(/^\s*\\ref\s*\{([^\}]+)\}\s*$/);
    if (!result) {return 0}
    var ref = result[1];
    if (this._labels[ref]) {
      this.setRef(element,ref);
    } else {
      if (!this._refs[ref]) {this._refs[ref] = []}
      this._refs[ref][this._refs[ref].length] = element;
    }
    return 1;
  },

  setRef: function (element,ref) {
    element.innerHTML = this.formatRef(this._labels[ref]);
    element.className = "label-ref";
  },

  fixRefs: function (label) {
    for (var i = 0; i < this._refs[label].length; i++)
      {this.setRef(this._refs[label][i],label)}
    delete this._refs[label];
  },

  badRefs: function () {
    for (var label in this._refs) {
      for (var i = 0; i < this._refs[label].length; i++) {
        var element = this._refs[label][i];
        element.className = "typeset";
        element.innerHTML = "<span class='error'>Reference '"+label+"' is undefined</span>";
      }
    }
  },

  makeDIV: function (element) {
    var div = document.createElement('div');
    div.className = 'jsMath_displayBox';
    div.innerHTML = '<div class="jsMath_number">' + this.nextNumber() + '</div>';
    element.parentNode.insertBefore(div,element);
    element.parentNode.removeChild(element);
    div.appendChild(element);
  },

  makeSPAN: function (element) {
    var span = document.createElement('span');
    span.className = 'jsMath_number';
    span.style.display = 'inline-block';
    span.innerHTML = jsMath.EqnNumber.nextNumber();
    element.parentNode.insertBefore(span,element);
  },

  ConvertMath: function (style,element,nocache) {
    var EqnNumber = jsMath.EqnNumber;
    if (EqnNumber.isRef(element)) return;
    EqnNumber._label = null; EqnNumber._nolabel = 0;
    this.ConvertMath_old(style,element,nocache);
    if (EqnNumber._label || (EqnNumber.autonumber && !EqnNumber._nolabel)) {
      if (element.tagName.toLowerCase() == 'div') {
        EqnNumber.makeDIV(element);
      } else if (element.parentNode.className == 'tex2math_div') {
        EqnNumber.makeSPAN(element);
      }
    }
  },

  ProcessComplete: function () {
    jsMath.EqnNumber.badRefs();
    this.ProcessComplete_old.apply(this,arguments);
  },

  Init: function () {
    jsMath.Setup.Styles(this.styles);
    jsMath.Translate.ConvertMath_old = jsMath.Translate.ConvertMath;
    jsMath.Translate.ConvertMath = this.ConvertMath;
    jsMath.Translate.ProcessComplete_old = jsMath.Translate.ProcessComplete;
    jsMath.Translate.ProcessComplete = this.ProcessComplete;
  },

  environments: {
    'equation*': 'Star',
    'eqnarray*': 'Star',
    'align*':    'Star',
    'multline*': 'Star',
    'gather*':   'Star',
    align:       ['StarExtension','AMSmath'],
    multline:    ['StarExtension','AMSmath'],
    gather:      ['StarExtension','AMSmath']
  },

  ResetStarEnvironments: function () {
    var Nenv = jsMath.EqnNumber.environments;
    var Penv = jsMath.Parser.prototype.environments;
    for (var name in Nenv) {
      if (name.match(/\*$/)) {Penv[name] = Nenv[name]}
    }
  }

};

if (jsMath.EqnNumber_old) {
  jsMath.Insert(jsMath.EqnNumber,jsMath.EqnNumber_old);
  delete jsMath.EqnNumber_old;
}

jsMath.Package(jsMath.Parser,{
  macros: {
    label:    'Label',
    nolabel:  'NoLabel',
    nonumber: 'NoLabel',
    ref:      'Ref'
  },

  environments: jsMath.EqnNumber.environments,

  Label: function (name) {
    var label = this.GetArgument(this.cmd+name); if (this.error) return;
    var EqnNumber = jsMath.EqnNumber;
    if (!EqnNumber._label) {
      if (!EqnNumber._labels[label]) {
        EqnNumber._label = label;
        EqnNumber._nolabel = 0;
      } else {
        this.Error("Label '"+label+"' is already defined");
      }
    } else {
      this.Error(this.cmd+name+' can only be used once in an equation');
    }
  },

  NoLabel: function (name) {
    var EqnNumber = jsMath.EqnNumber;
    EqnNumber._label = null; EqnNumber._nolabel = 1;
  },

  Ref: function (name) {
    this.Error(this.cmd+name+' must be used by itself');
  },

  Star: function (name) {
    this.NoLabel();
    var cmd = this.environments[name.substr(0,name.length-1)];
    if (typeof(cmd) === 'string') {cmd = [cmd]}
    this[cmd[0]](name,cmd.slice(1));
  },

  StarExtension: function (name,data) {
    try {this.Extension(name,data)} catch (e) {}
    jsMath.Synchronize(jsMath.EqnNumber.ResetStarEnvironments);
    throw "restart";
  }

});

jsMath.EqnNumber.Init();

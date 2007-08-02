/*
 *  extensions/newcommand.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file implements the \newcommand and \def macros.  It will be
 *  loaded automatically when needed, or can be loaded by
 *
 *    jsMath.Extension.Require('newcommand');
 *
 *  ---------------------------------------------------------------------
 *
 *  Copyright 2005-2006 by Davide P. Cervone
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

jsMath.Package(jsMath.Parser,{

  macros: {
    newcommand: 'NewCommand',
    newenvironment: 'NewEnvironment',
    def: 'MacroDef'
  },

  /*
   *  Implement \newcommand{\name}[n]{...}
   */
  NewCommand: function (name) {
    var cs = this.trimSpaces(this.GetArgument(this.cmd+name)); if (this.error) return;
    var n = this.trimSpaces(this.GetBrackets(this.cmd+name)); if (this.error) return;
    var def = this.GetArgument(this.cmd+name); if (this.error) return;
    if (n == '') {n = null}
    if (cs.charAt(0) == this.cmd) {cs = cs.substr(1)}
    if (!cs.match(/^(.|[a-z]+)$/)) {this.Error("Illegal control sequence name for "+this.cmd+name); return}
    if (n != null && !n.match(/^[0-9]+$/)) {this.Error("Illegal number of parameters specified in "+this.cmd+name); return}
    jsMath.Parser.prototype.macros[cs] = ['Macro',def,n];
  },

  /*
   *  Implement \newenvironment{name}[n]{begincmd}{endcmd}
   */
  NewEnvironment: function (name) {
    var env = this.trimSpaces(this.GetArgument(this.cmd+name)); if (this.error) return;
    var n = this.trimSpaces(this.GetBrackets(this.cmd+name)); if (this.error) return;
    var bdef = this.GetArgument(this.cmd+name); if (this.error) return;
    var edef = this.GetArgument(this.cmd+name); if (this.error) return;
    if (n == '') {n = null}
    if (n != null && !n.match(/^[0-9]+$/)) {this.Error("Illegal number of parameters specified in "+this.cmd+name); return}
    jsMath.Parser.prototype.environments[env] = ['Environment',bdef,edef,n];
  },

  /*
   *  Implement \def command
   */
  MacroDef: function (name) {
    var cs = this.GetCSname(this.cmd+name); if (this.error) return;
    var params = this.GetTemplate(this.cmd+name); if (this.error) return;
    var def = this.GetArgument(this.cmd+name); if (this.error) return;
    if (typeof(params) == 'number') {
      jsMath.Parser.prototype.macros[cs] = ['Macro',def,params];
    } else {
      jsMath.Parser.prototype.macros[cs] = ['MacroWithTemplate',def,params[0],params[1]];
    }
  },

  /*
   *  Get a CS name or give an error
   */
  GetCSname: function (cmd) {
    var c = this.GetNext();
    if (c != this.cmd) {this.Error(cmd+" must be followed by a control sequence"); return null}
    var cs = this.trimSpaces(this.GetArgument(cmd)); if (this.error) {return null};
    return cs.substr(1);
  },

  /*
   *  Get a \def parameter template
   */
  GetTemplate: function (cmd) {
    var c; var params = []; var n = 0;
    c = this.GetNext(); var i = this.i;
    while (this.i < this.string.length) {
      c = this.GetNext();
      if (c == '#') {
        if (i != this.i) {params[n] = this.string.substr(i,this.i-i)}
        c = this.string.charAt(++this.i);
        if (!c.match(/[1-9]/)) {this.Error("Illegal use of # in "+cmd); return null}
        if (1*c != ++n) {this.Error("Parameters must be numbered sequentially"); return null}
        i = this.i+1;
      } else if (c == '{') {
        if (i != this.i) {params[n] = this.string.substr(i,this.i-i)}
        if (params.length > 0) {return [n,params]} else {return n}
      }
      this.i++;
    }
    this.Error("Missing replacement string for definition of "+cmd);
    return null;
  },

  /*
   *  Process a macro with a parameter template
   */
  MacroWithTemplate: function (name,data) {
    var text = data[0];
    var n = data[1]; var params = data[2];
    if (n) {
      var args = []; var c = this.GetNext();
      if (params[0] && !this.MatchParam(params[0]))
        {this.Error("Use of "+this.cmd+name+" doesn't match its definition"); return}
      for (var i = 0; i < n; i++) {
        args[args.length] = this.GetParameter(this.cmd+name,params[i+1]);
        if (this.error) return;
      }
      text = this.SubstituteArgs(args,text);
    }
    this.string = this.AddArgs(text,this.string.slice(this.i));
    this.i = 0;
  },

  /*
   *  Process a user-defined environment
   */
  Environment: function (name,data) {
    var bdef = data[0]; var edef = data[1]; var n = data[2];
    if (n) {
      var args = [];
      for (var i = 0; i < n; i++) {
        args[args.length] = this.GetArgument(this.cmd+"begin{"+name+"}"); if (this.error) return;
      }
      bdef = this.SubstituteArgs(args,bdef);
    }
    var text = this.GetEnd(name); if (this.error) return;
    text = this.AddArgs(this.AddArgs(bdef,text),edef);
    this.string = this.AddArgs(text,this.string.slice(this.i));
    this.i = 0;
  },

  /*
   *  Find a single parameter delimited by a trailing template
   */
  GetParameter: function (name,param) {
    if (param == null) {return this.GetArgument(name)}
    var i = this.i; var j = 0; var hasBraces = 0;
    while (this.i < this.string.length) {
      if (this.string.charAt(this.i) == '{') {
        if (this.i == i) {hasBraces = 1}
        this.GetArgument(name); j = this.i - i;
      } else if (this.MatchParam(param)) {
        if (hasBraces) {i++; j -= 2}
        return this.string.substr(i,j);
      } else {
        this.i++; j++; hasBraces = 0;
      }
    }
    this.Error("Runaway argument for "+name+"?");
    return null;
  },

  /*
   *  Check if a template is at the current location.
   *  (The match must be exact, with no spacing differences.  TeX is
   *   a little more forgiving about spaces after macro names)
   */
  MatchParam: function (param) {
    if (this.string.substr(this.i,param.length) != param) {return 0}
    this.i += param.length;
    return 1;
  }

});

/*
 *  Define a jsMath.Environment() command similar to the
 *  jsMath.Macro() command.
 *
 *  Usage:  jsMath.Environment(name,begin,end[,n])
 *
 *  where "name" is the name of the environment, "begin" is the
 *  text that replaces the \begin{name} and "end" is the text that
 *  replaces the \end{name}.  If "n" is provided, it is the number
 *  of parameters that the \begin{name} accepts, and these are
 *  used to replace #1, #2, etc within the "begin" text.
 */

jsMath.Add(jsMath,{
  Environment: function (name) {
    var environments = jsMath.Parser.prototype.environments;
    environments[name] = ['Environment'];
    for (var i = 1; i < arguments.length; i++)
      {environments[name][environments[name].length] = arguments[i]}
  }
});

jsMath.Package(jsMath.Parser,{

  macros: {font: 'Font'},
  fontCS: {},

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
   *  Handle the \font command
   */
  Font: function (name) {
    var cs = this.GetCSname(this.cmd+name); if (this.error) return;
    while (this.nextIsSpace()) {this.i++}
    if (this.string.charAt(this.i++) == '=') {
      while (this.nextIsSpace()) {this.i++}
      var font = this.string.slice(this.i).match(/^[a-z]+[0-9]+/i);
      if (font) {
        this.i += (new String(font)).length;
        if (jsMath.TeX.famName[font] != null) {
          this.macros[cs] = ['HandleFont',jsMath.TeX.famName[font]];
        } else {
          this.macros[cs] = ['Extension',jsMath.Font.URL(font),"fontCS"];
          this.fontCS[cs] = 1; // so Extension has something to delete
        }
      } else {this.Error("Missing font name")}
    } else {this.Error("Missing font definition")}
  }

});


function checkForGearsInstalled() {
    if (window.google && window.google.gears) {
        var controls = document.getElementById("controls");
        var newSpan = document.createElement("span");
        newSpan.setAttribute("class","vbar");
        var newLink = document.createElement("a");
        newLink.setAttribute("class","usercontrol");
        newLink.setAttribute("onClick","createShortcut()");
        newLink.setAttribute("title","Create application shortcut with gears");
        var newText = document.createTextNode("Create Shortcut");
        newLink.appendChild(newText)
        controls.appendChild(newSpan)
        controls.appendChild(newLink)
    }
}

function createShortcut() {
   var desktop = google.gears.factory.create('beta.desktop');
   desktop.createShortcut('SAGE',
       '.',
       {'128x128': '../../images/icon128x128.png',
           '48x48': '../../images/icon48x48.png',
           '32x32': '../../images/icon32x32.png',
           '16x16': '../../images/icon16x16.png'},
           'SAGE Notebook');
}

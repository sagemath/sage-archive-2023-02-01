package org.sagemath.sage3d;

import java.awt.event.WindowEvent;
import java.awt.event.WindowAdapter;
import java.io.File;
import java.net.URL;
import java.net.MalformedURLException;
import javax.swing.JFrame;


public class ObjectViewerApp extends ObjectViewerApplet {

  public static void main(String[] args) {
    String name = null, id = "", url = args[0];

    if (args.length > 1) {
      name = args[1];
    }
    else {
      int ix = url.lastIndexOf('/');
      if (ix == -1) name = url;
      else name = url.substring(ix);
    }

    if (args.length > 2) {
      id = args[2];
    }
    new ObjectViewerApp().showView(url, id, name);
  }

  public URL getBase() {
    try {
      return new File(".").toURL();
    }
    catch (MalformedURLException ex) {
      System.err.println(ex);
      return null;
    }
  }

  public JFrame getWindow(String id, String name) {
    JFrame frame = super.getWindow(id, name);
    frame.addWindowListener( new WindowAdapter() {
      public void windowClosing(WindowEvent e) {
        System.exit(0);
      }
    });
    return frame;
  }

}
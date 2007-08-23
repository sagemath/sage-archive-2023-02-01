package org.sagemath.sage3d;

import java.applet.Applet;
import javax.swing.JApplet;
import javax.swing.JPanel;
import java.awt.BorderLayout;
import javax.swing.JComponent;

import java.util.HashMap;

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.net.*;


import java.util.Vector;



//import com.sun.j3d.loaders.objectfile.*;


//import java.net.MalformedURLException;

//import java.io.FileNotFoundException;
import javax.swing.*;

//import com.sun.j3d.utils.applet.MainFrame;

public class ObjectViewerApplet extends JApplet {

  protected HashMap windows = new HashMap();
  protected ObjectViewer curView;

  public void start() {
    setLayout(new BorderLayout());
    String inline = getParameter("inline");
    if (inline != null) showView(inline, null);
    String popup = getParameter("popup");
    if (popup != null) showView(popup, "popup");
  }

  public void showView(String url, String id) {
    showView(url, id, url);
  }

  public void showView(URL url, String id) {
    showView(url, id, url.toString());
  }

  public void showView(String url, String id, String name) {
    try {
      showView(new URL(getDocumentBase(), url), id, name);
    }
    catch (MalformedURLException ex) {
      System.out.println(ex);
      showView((URL)null, null, null);
    }
  }

  public void showView(URL url, String id, String name) {
    ObjectViewer view = new ObjectViewer(url);
    if (id == null) {
      getContentPane().add(view);
      validate();
      curView = view;
    }
    else {
      JFrame frame = (JFrame)windows.get(id);
      if (frame == null) {
        frame = new JFrame(name);
        frame.getContentPane().setLayout(new BorderLayout());
        //windows.put(id, frame); // does holding on to it cause things to crash?
        frame.setSize(300,300);
      }
      else {
        frame.getContentPane().removeAll();
      }
      frame.getContentPane().add(view);
      frame.setVisible(true);
    }
  }

}

package org.sagemath.sage3d;

import javax.swing.JApplet;
import java.awt.BorderLayout;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

import java.io.*;
import java.net.*;

import java.util.Vector;
import java.util.HashMap;


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
      showView(new URL(getBase(), url), id, name);
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
      JFrame frame = getWindow(id, name);
      frame.getContentPane().add(view);
      frame.setVisible(true);
    }
  }

  public URL getBase() {
    return getDocumentBase();
  }

  public JFrame getWindow(String id, String name) {
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
    return frame;
  }

}

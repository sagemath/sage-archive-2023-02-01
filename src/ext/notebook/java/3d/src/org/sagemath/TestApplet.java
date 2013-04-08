package org.sagemath;

import java.applet.*;
import java.awt.*;

import javax.swing.JFrame;
import java.util.HashMap;

public class TestApplet extends Applet {

  protected HashMap windows = new HashMap();

  public TestApplet() {
    setLayout(new BorderLayout());
    add(new Label("Hello SAGE"));
  }

  public void showView(String url, String id) {
    showView(url, id, url);
  }

  public void showView(String url, String id, String name) {
    JFrame frame = (JFrame)windows.get(id);
    if (frame == null) {
      frame = new JFrame(id);
      frame.setLayout(new BorderLayout());
      windows.put(id, frame);
    }
    frame.getContentPane().add(new Label(url));
    frame.setSize(200,50);
    frame.setVisible(true);
  }

}
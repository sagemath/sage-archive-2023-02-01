import javax.swing.JPanel;
import java.awt.*;

import java.util.HashMap;
import java.util.Enumeration;
import java.io.*;
import java.net.URL;

import com.sun.j3d.utils.universe.*;
import com.sun.j3d.utils.geometry.*;
import javax.media.j3d.*;
import javax.vecmath.*;

import com.sun.j3d.loaders.*;
import com.sun.j3d.loaders.objectfile.ObjectFile;

import com.sun.j3d.utils.behaviors.keyboard.*;
import com.sun.j3d.utils.behaviors.mouse.*;



public class ObjectViewer extends JPanel {

//  protected SimpleUniverse univ;
  protected BranchGroup scene;
  protected TransformGroup transform;
  protected TransformGroup viewStart;

  protected Canvas3D canvas;



  public ObjectViewer(URL url) {
    setLayout(new BorderLayout());
    Canvas3D canvas3D = new Canvas3D(SimpleUniverse.getPreferredConfiguration());
    add("Center", canvas3D);

    BoundingSphere bounds = new BoundingSphere(new Point3d(), 1000);
    BranchGroup root = new BranchGroup();
    BranchGroup scene = createSceneGraph(url);
    scene.setBoundsAutoCompute(true);
    System.out.println(scene.getBounds());
    BoundingSphere sceneBounds = new BoundingSphere(scene.getBounds());

    SimpleUniverse univ = new SimpleUniverse(canvas3D);
    ViewingPlatform view = univ.getViewingPlatform();
    view.setNominalViewingTransform();

    Transform3D t = new Transform3D();
    TransformGroup viewTransform = view.getViewPlatformTransform();

    t.set(new Vector3d(0,0,5*sceneBounds.getRadius()));
    viewTransform.setTransform(t);

    BranchGroup lights = new BranchGroup();
    Light light = new AmbientLight();
    light.setInfluencingBounds(bounds);
    lights.addChild(light);
    light = new DirectionalLight();
    light.setInfluencingBounds(bounds);
    lights.addChild(light);
    root.addChild(lights);

    TransformGroup tg = new TransformGroup();
    tg.setCapability(TransformGroup.ALLOW_TRANSFORM_READ);
    tg.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
    tg.addChild(scene);

    root.addChild(tg);

    MouseRotate mouse = new MouseRotate();
    mouse.setTransformGroup(tg);
    mouse.setSchedulingBounds(bounds);
    root.addChild(mouse);

    MouseZoom mousezoom = new MouseZoom();
    mousezoom.setTransformGroup(tg);
    mousezoom.setSchedulingBounds(bounds);
    root.addChild(mousezoom);
    root.compile();

    univ.addBranchGraph(root);

  }


  public BranchGroup createSceneGraph(URL url) {

    try {
      Loader loader = new com.sun.j3d.loaders.objectfile.ObjectFile();
      Scene scene  = loader.load(url);
      BranchGroup bg = scene.getSceneGroup();
      System.out.println(bg);
      TransformGroup[] views = scene.getViewGroups();
      if (views != null) {
        for(int i=0; i<views.length; i++) {
          System.out.print(views[i]);
        }
        if (views.length > 0) viewStart = views[0];
      }
      return bg;
    }
    catch (Exception ex) {
      //in case there was a problem, print the stack out
      //ex.printStackTrace();
      System.out.println(ex);
      BranchGroup bg = new BranchGroup();
      bg.addChild(new ColorCube());
      System.out.println(bg);
      return bg;
    }

  }

}

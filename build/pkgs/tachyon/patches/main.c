/*
 * main.c - This file contains the main program and driver for the raytracer.
 *
 *  $Id: main.c,v 1.76 2010/01/18 19:36:34 johns Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tachyon.h"    /* The Tachyon ray tracing library API */
#include "getargs.h"    /* command line argument/option parsing */
#include "parse.h"      /* Support for my own scene file format */
#include "nffparse.h"   /* Support for NFF files, as in SPD */
#include "ac3dparse.h"  /* Support for AC3D files */
#include "mgfparse.h"   /* Support for MGF files */

#ifdef USEOPENGL
#include "glwin.h"      /* OpenGL run-time display code */
#endif

#ifdef USERTVI
#include "rtvi_iface.h" /* Synergy RTVI/ETVI Attached Framebuffers */
#endif

#ifdef USESPACEBALL
#include "spaceball.h"  /* Spaceball fly-through code */
#endif

typedef struct {
  float x;
  float y;
  float z;
} floatvec;


typedef struct {
  int xsize, ysize;

#if defined(USEOPENGL) || defined(USERTVI)
    unsigned char * img;
#endif
#ifdef USEOPENGL
    void * glwin;
#endif
#ifdef USERTVI
    void * rtviwin;
#endif
} dispHandle;


static void my_ui_message(int a, char * msg) {
  printf("%s\n", msg);
}

static void my_ui_progress(int percent) {
  printf("\rRendering Progress:       %3d%% complete            \r", percent);
  fflush(stdout);
}


/*
 * routines for managing runtime display of ray traced scene
 */
static dispHandle * tachyon_display_create(SceneHandle scene) {
  dispHandle * dh;

  dh = (dispHandle *) malloc(sizeof(dispHandle));

  if (dh != NULL) {
    memset(dh, 0, sizeof(dispHandle));

    rt_get_resolution(scene, &dh->xsize, &dh->ysize);

#if defined(USEOPENGL) || defined(USERTVI)
      dh->img = malloc((dh->xsize)*(dh->ysize)*3);
      if (dh->img != NULL) {

#if defined(USEOPENGL)
        dh->glwin = glwin_create("Tachyon Parallel/Multiprocessor Ray Tracer", dh->xsize, dh->ysize);
#elif defined(USERTVI)
        dh->rtviwin = rt_rtvi_init(dh->xsize, dh->ysize);
#endif

        rt_rawimage_rgb24(scene, dh->img);
      }
      else {
        printf("Couldn't allocate image buffer for framebuffer display!!\n");
        free(dh);
        return NULL;
      }
#endif
  }

  return dh;
}

static void tachyon_display_draw(dispHandle *dh) {
#if defined(USEOPENGL)
      if (dh->img != NULL) {
        glwin_handle_events(dh->glwin);
        glwin_draw_image(dh->glwin, dh->xsize, dh->ysize, dh->img);
      }
#elif defined(USERTVI)
      if (dh->img != NULL)
        rt_rtvi_displayimage(dh->img, rtviwin);
#endif
}

static void tachyon_display_delete(dispHandle *dh) {
#if defined(USEOPENGL) || defined(USERTVI)
    if (dh->img != NULL) {
#if defined(USEOPENGL)
      glwin_destroy(dh->glwin);
#endif
      free(dh->img);
    }
#endif
}


/*
 * main loop for creating animations by flying using a spaceball
 * or other 3-D input mechanism.
 */
static int fly_scene(argoptions opt, SceneHandle scene, int node) {
  dispHandle * dh = NULL;
  int done = 0;
  int frameno = 0;
  float fps;
  rt_timerhandle fpstimer;
  rt_timerhandle animationtimer;
  char outfilename[1];

#if defined(USESPACEBALL)
  sbHandle * bh = NULL;
#endif

  if (node == 0)
    dh = tachyon_display_create(scene);

  rt_set_ui_message(NULL);
  rt_set_ui_progress(NULL);

  if (node == 0)
    printf("Interactive Camera Flight\n");

  outfilename[0] = '\0';
  rt_outputfile(scene, outfilename);

  fpstimer=rt_timer_create();
  animationtimer=rt_timer_create();

#if defined(USESPACEBALL)
  if (node == 0) {
#if 1
    bh = tachyon_init_spaceball(scene, opt.spaceball);
#else
    if (rt_numnodes() < 2) {
      bh = tachyon_init_spaceball(scene, opt.spaceball);
    } else {
      printf("WARNING: Spaceball mode disabled when running with distributed memory");
    }
#endif
  }
#endif

  rt_timer_start(animationtimer);
  while (!done) {
    if (frameno != 0) {
      rt_timer_stop(fpstimer);
      fps = 1.0f / rt_timer_time(fpstimer);
    } else {
      fps = 0.0;
    }

    rt_timer_start(fpstimer);
    if (node == 0) {
      printf("\rRendering Frame: %9d   %10.4f FPS       ", frameno, fps);
      fflush(stdout);
    }

#if defined(USESPACEBALL)
    if (bh != NULL)
      done = tachyon_spaceball_update(bh, scene);
#endif

    rt_renderscene(scene);

    if (dh != NULL)
      tachyon_display_draw(dh);

    frameno++;
  }

  rt_timer_stop(animationtimer);
  fps = frameno / rt_timer_time(animationtimer);

  if (node == 0) {
    printf("\rCompleted animation of %d frames                            \n", frameno);
    printf("Animation Time: %10.4f seconds  (Averaged %7.4f FPS)\n",
         rt_timer_time(animationtimer), fps);
  }
  rt_timer_destroy(fpstimer);

  if (node == 0) {
    printf("\nFinished Running Camera.\n");

    if (dh !=NULL)
      tachyon_display_delete(dh);
  }

  rt_deletescene(scene); /* free the scene */
  rt_finalize(); /* close down the rendering library and MPI */

  return 0;
}



/*
 * main loop for creating animations by playing recorded camera fly-throughs
 */
static int animate_scene(argoptions opt, SceneHandle scene, int node) {
  char outfilename[1000];
  FILE * camfp;
  dispHandle * dh = NULL;

  if (node == 0)
    dh = tachyon_display_create(scene);

  /* if we have a camera file, then animate.. */
  if ((camfp = fopen(opt.camfilename, "r")) != NULL) {
    floatvec cv, cu, cc;
    apivector cmv, cmu, cmc;
    int frameno = 0;
    float fps;
    rt_timerhandle fpstimer;
    rt_timerhandle animationtimer;

    rt_set_ui_message(NULL);
    rt_set_ui_progress(NULL);

    if (node == 0)
      printf("Running Camera File: %s\n", opt.camfilename);

    fpstimer=rt_timer_create();
    animationtimer=rt_timer_create();

    rt_timer_start(animationtimer);

    while (!feof(camfp)) {
      fscanf(camfp, "%f %f %f  %f %f %f  %f %f %f",
        &cv.x, &cv.y, &cv.z, &cu.x, &cu.y, &cu.z, &cc.x, &cc.y, &cc.z);

      cmv.x = cv.x; cmv.y = cv.y; cmv.z = cv.z;
      cmu.x = cu.x; cmu.y = cu.y; cmu.z = cu.z;
      cmc.x = cc.x; cmc.y = cc.y; cmc.z = cc.z;

      if (frameno != 0) {
        rt_timer_stop(fpstimer);
        fps = 1.0f / rt_timer_time(fpstimer);
      } else {
        fps = 0.0;
      }

      rt_timer_start(fpstimer);
      outfilename[0] = '\0';
      if (opt.nosave == 1) {
        if (node == 0) {
          printf("\rRendering Frame: %9d   %10.4f FPS       ", frameno, fps);
          fflush(stdout);
        }
      }
      else {
        sprintf(outfilename, opt.outfilename, frameno);
        if (node == 0) {
          printf("\rRendering Frame to %s   (%10.4f FPS)       ", outfilename, fps);
          fflush(stdout);
        }
      }

      rt_outputfile(scene, outfilename);
      rt_camera_position(scene, cmc, cmv, cmu);

      rt_renderscene(scene);

      if (dh != NULL)
        tachyon_display_draw(dh);

      frameno++;
    }
    rt_timer_stop(animationtimer);
    fps = frameno / rt_timer_time(animationtimer);
    if (node == 0) {
      printf("\rCompleted animation of %d frames                            \n", frameno);
      printf("Animation Time: %10.4f seconds  (Averaged %7.4f FPS)\n",
           rt_timer_time(animationtimer), fps);
    }
    rt_timer_destroy(fpstimer);
    fclose(camfp);
  } else {
    if (node == 0) {
      printf("Couldn't open camera file: %s\n", opt.camfilename);
      printf("Aborting render.\n");
    }
    rt_deletescene(scene); /* free the scene */
    rt_finalize(); /* close down the rendering library and MPI */
    return -1;
  }

  if (node == 0) {
    printf("\nFinished Running Camera.\n");

    if (dh !=NULL)
      tachyon_display_delete(dh);
  }

  rt_deletescene(scene); /* free the scene */
  rt_finalize(); /* close down the rendering library and MPI */

  return 0;
}




#ifdef VXWORKS
int ray(int argc, char **argv) {
#else
int main(int argc, char **argv) {
#endif
  SceneHandle scene;
  unsigned int rc;
  argoptions opt;
  char * filename;
  int node, fileindex;
  rt_timerhandle parsetimer;
  size_t len;

  node = rt_initialize(&argc, &argv);

  rt_set_ui_message(my_ui_message);
  rt_set_ui_progress(my_ui_progress);

  if (node == 0) {
    printf("Tachyon Parallel/Multiprocessor Ray Tracer   Version %s   \n",
           TACHYON_VERSION_STRING);
    printf("Copyright 1994-2010,    John E. Stone <john.stone@gmail.com> \n");
    printf("------------------------------------------------------------ \n");
  }

  if ((rc = getargs(argc, argv, &opt, node)) != 0) {
    rt_finalize();
    exit(rc);
  }

  if (opt.numfiles > 1) {
    printf("Rendering %d scene files.\n", opt.numfiles);
  }

  for (fileindex=0; fileindex<opt.numfiles; fileindex++) {
    scene = rt_newscene();

    /* process command line overrides */
    presceneoptions(&opt, scene);

    filename = opt.filenames[fileindex];

    if (opt.numfiles > 1) {
      printf("\nRendering scene file %d of %d, %s\n", fileindex+1, opt.numfiles, filename);
    }

    parsetimer=rt_timer_create();
    rt_timer_start(parsetimer);

    len = strlen(filename);

    if (len > 4 && (!strcmp(filename+len-4, ".nff") ||
                    !strcmp(filename+len-4, ".NFF"))) {
      rc = ParseNFF(filename, scene); /* must be an NFF file */
    }
    else if (len > 3 && (!strcmp(filename+len-3, ".ac") ||
                         !strcmp(filename+len-3, ".AC"))) {
      rc = ParseAC3D(filename, scene); /* Must be an AC3D file */
    }
#ifdef USELIBMGF
    else if (len > 4 && (!strcmp(filename+len-4, ".mgf") ||
                         !strcmp(filename+len-4, ".MGF"))) {
      rc = ParseMGF(filename, scene, 1); /* Must be an MGF file */
    }
#endif
    else {
      rc = readmodel(filename, scene); /* Assume its a Tachyon scene file */
    }

    rt_timer_stop(parsetimer);
    if (rc == PARSENOERR && node == 0)
      printf("Scene Parsing Time: %10.4f seconds\n", rt_timer_time(parsetimer));
    rt_timer_destroy(parsetimer);

    if (rc != PARSENOERR && node == 0) {
      switch(rc) {
        case PARSEBADFILE:
          printf("Parser failed due to nonexistent input file: %s\n", filename);
          break;
        case PARSEBADSUBFILE:
          printf("Parser failed due to nonexistent included file.\n");
          break;
        case PARSEBADSYNTAX:
          printf("Parser failed due to an input file syntax error.\n");
          break;
        case PARSEEOF:
          printf("Parser unexpectedly hit an end of file.\n");
          break;
        case PARSEALLOCERR:
          printf("Parser ran out of memory.\n");
          break;
      }
      if (fileindex+1 < opt.numfiles)
        printf("Aborting render, continuing with next scene file...\n");
      else
        printf("Aborting render.\n");

      rt_deletescene(scene); /* free the scene */
      continue;              /* process the next scene */
    }

    /* process command line overrides */
    postsceneoptions(&opt, scene);

    /* choose which rendering mode to use */
    if (opt.usecamfile == 1) {
      return animate_scene(opt, scene, node); /* fly using prerecorded data */
    }
    else if (strlen(opt.spaceball) > 0) {
      return fly_scene(opt, scene, node);     /* fly with spaceball etc */
    }
    else {
      if (opt.numfiles > 1 && opt.nosave != 1) {
        char multioutfilename[FILENAME_MAX];
        sprintf(multioutfilename, opt.outfilename, fileindex);
        rt_outputfile(scene, multioutfilename);
      }

      rt_renderscene(scene); /* Render a single frame */
    }

    rt_deletescene(scene);   /* free the scene, get ready for next one */
  }

  rt_finalize();             /* close down the rendering library and MPI */
  freeoptions(&opt);         /* free parsed command line option data */

  return 0;
}



//
//  PreferencePanelController.h
//
//  Created by Ivan Andrus on 26/6/10.
//  Copyright 2010 __MyCompanyName__. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface PreferencePanelController : NSObject {
    IBOutlet id appController;
    IBOutlet id prefWindow;

    IBOutlet id showInDock;
    IBOutlet id TerminalApplescript;
    IBOutlet id TerminalEmulator;
    IBOutlet id SessionType;
    IBOutlet id DefaultArgs;
}

- (IBAction)apply:(id)sender;
- (IBAction)resetTerminalApplescript:(id)sender;
- (IBAction)addToPATH:(id)sender;
- (void)windowDidBecomeKey:(NSNotification *)aNotification;
- (void)comboBoxWillDismiss:(NSNotification *)notification;
- (void)controlTextDidEndEditing:(NSNotification *)aNotification;
- (void)updateForComboBoxChanges;

@end

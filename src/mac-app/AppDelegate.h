//
//  AppDelegate.h
//
//  Created by Ivan Andrus on 26/6/10.
//  Copyright 2010 __MyCompanyName__. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import "AppController.h"

// http://stackoverflow.com/questions/1496788/building-for-10-5-in-xcode-3-2-on-snow-leopard-error

// this is a 10.6 only protocol, but we build on 10.4 and 10.5 which
// makes the #if more complicated than the webpage, so we don't worrry
// about it.
//@interface AppDelegate : NSObject <NSApplicationDelegate>
@interface AppDelegate : NSObject
{
    IBOutlet AppController* appController;

}

+ (void)initialize;
- (void)applicationDidFinishLaunching:(NSNotification *)aNotification;
- (void)applicationWillFinishLaunching:(NSNotification *)aNotification;
- (BOOL)applicationShouldOpenUntitledFile:(NSApplication *)sender;

// Opening files
- (BOOL)application:(NSApplication * )theApplication openFile: (NSString * )filename;
- (void)getUrl:(NSAppleEventDescriptor *)event withReplyEvent:(NSAppleEventDescriptor *)replyEvent;


-(IBAction)openDocumentWithDialogBox:(id)sender;

@end

//
//  InputPanelController.m
//
//  Created by Ivan Andrus on 7/9/10.
//  Copyright 2010 __MyCompanyName__. All rights reserved.
//

#import "InputPanelController.h"
#import "AppController.h"

@implementation InputPanelController

- (void)runCommand:(NSString*)command withPrompt:(NSString*)prompt withArguments:(NSString*)defArgs editingCommand:(BOOL)editable{

    // If it wasn't closed, give the user a bit of a flicker so they know it's a different prompt
    [self close];

    [label setStringValue:prompt];
    if (editable) {
        [textField setStringValue: (defArgs == nil) ? command : [command stringByAppendingFormat:@" %@", defArgs]];
        commandPrefix = @"";
    } else {
        [textField setStringValue: (defArgs == nil ) ? @"" : defArgs];
        commandPrefix = [command retain];
    }
    [NSApp activateIgnoringOtherApps:YES];
    [window makeKeyAndOrderFront:self];
    [self showWindow:self];
}

- (IBAction)accept:(id)sender{
    [self close];
    NSString * command = [NSString stringWithFormat:@"%@ %@", commandPrefix, [textField stringValue]];
    [appController terminalRun:command];
    [commandPrefix release];
}


@end

//
//  InputPanelController.h
//
//  Created by Ivan Andrus on 7/9/10.
//  Copyright 2010 __MyCompanyName__. All rights reserved.
//

#import <Cocoa/Cocoa.h>

@interface InputPanelController : NSWindowController {
    IBOutlet id label;
    IBOutlet id textField;
    IBOutlet id window;
    IBOutlet id appController;

    NSString *commandPrefix;
}

- (void)runCommand:(NSString*)command withPrompt:(NSString*)prompt withArguments:(NSString*)defArgs editingCommand:(BOOL)editable;
- (IBAction)accept:(id)sender;

@end

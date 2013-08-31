//
//  MyDocument.h
//  Sage
//
//  Created by Ivan Andrus on 26/6/10.
//  Copyright 2010 __MyCompanyName__. All rights reserved.
//


#import <Cocoa/Cocoa.h>
#import <WebKit/WebView.h>

@interface MyDocument : NSDocument
{
    IBOutlet id urlString;
    IBOutlet WebView *webView;
    IBOutlet id progressIndicator;
}

- (IBAction)connectURL:(id)sender;
- (id)webView;
- (void)webViewShow:(WebView *)sender;
- (WebView *)webView:(WebView *)sender createWebViewWithRequest:(NSURLRequest *)request;
- (IBAction)browseURL:(NSString*)theURL;
- (IBAction)printDocument:(id)sender;


@end

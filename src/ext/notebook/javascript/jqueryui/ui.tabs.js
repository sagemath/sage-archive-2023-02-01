/*
 * Tabs 3 - New Wave Tabs
 *
 * Copyright (c) 2007 Klaus Hartl (stilbuero.de)
 * Dual licensed under the MIT (MIT-LICENSE.txt)
 * and GPL (GPL-LICENSE.txt) licenses.
 *
 * http://docs.jquery.com/UI/Tabs
 */

(function($) {

    // if the UI scope is not availalable, add it
    $.ui = $.ui || {};

    // tabs API methods
    $.fn.tabs = function() {
        var method = typeof arguments[0] == 'string' && arguments[0];
        var args = method && Array.prototype.slice.call(arguments, 1) || arguments;

        return this.each(function() {
            if (method) {
                var tabs = $.data(this, 'ui-tabs');
                tabs[method].apply(tabs, args);
            } else
                new $.ui.tabs(this, args[0] || {});
        });
    };

    // tabs class
    $.ui.tabs = function(el, options) {
        var self = this;

        this.element = el;

        this.options = $.extend({

            // basic setup
            selected: 0,
            unselect: options.selected === null,
            event: 'click',
            disabled: [],
            cookie: null, // pass options object as expected by cookie plugin: { expires: 7, path: '/', domain: 'jquery.com', secure: true }
            // TODO bookmarkable: $.ajaxHistory ? true : false,

            // Ajax
            spinner: 'Loading&#8230;',
            cache: false,
            idPrefix: 'ui-tabs-',
            ajaxOptions: {},

            // animations
            fx: null, /* e.g. { height: 'toggle', opacity: 'toggle', duration: 200 } */

            // templates
            tabTemplate: '<li><a href="#{href}"><span>#{label}</span></a></li>',
            panelTemplate: '<div></div>',

            // CSS classes
            navClass: 'ui-tabs-nav',
            selectedClass: 'ui-tabs-selected',
            unselectClass: 'ui-tabs-unselect',
            disabledClass: 'ui-tabs-disabled',
            panelClass: 'ui-tabs-panel',
            hideClass: 'ui-tabs-hide',
            loadingClass: 'ui-tabs-loading'

        }, options);

        this.options.event += '.ui-tabs'; // namespace event
        this.options.cookie = $.cookie && $.cookie.constructor == Function && this.options.cookie;

        $(el).bind('setData.ui-tabs', function(event, key, value) {
            self.options[key] = value;
            this.tabify();
        }).bind('getData.ui-tabs', function(event, key) {
            return self.options[key];
        });

        // save instance for later
        $.data(el, 'ui-tabs', this);

        // create tabs
        this.tabify(true);
    };

    // instance methods
    $.extend($.ui.tabs.prototype, {
        tabId: function(a) {
            return a.title && a.title.replace(/\s/g, '_').replace(/[^A-Za-z0-9\-_:\.]/g, '')
                || this.options.idPrefix + $.data(a);
        },
        ui: function(tab, panel) {
            return {
                instance: this,
                options: this.options,
                tab: tab,
                panel: panel
            };
        },
        tabify: function(init) {

            this.$lis = $('li:has(a[href])', this.element);
            this.$tabs = this.$lis.map(function() { return $('a', this)[0]; });
            this.$panels = $([]);

            var self = this, o = this.options;

            this.$tabs.each(function(i, a) {
                // inline tab
                if (a.hash && a.hash.replace('#', '')) // Safari 2 reports '#' for an empty hash
                    self.$panels = self.$panels.add(a.hash);
                // remote tab
                else if ($(a).attr('href') != '#') { // prevent loading the page itself if href is just "#"
                    $.data(a, 'href.ui-tabs', a.href); // required for restore on destroy
                    $.data(a, 'load.ui-tabs', a.href); // mutable
                    var id = self.tabId(a);
                    a.href = '#' + id;
                    var $panel = $('#' + id);
                    if (!$panel.length) {
                        $panel = $(o.panelTemplate).attr('id', id).addClass(o.panelClass)
                            .insertAfter( self.$panels[i - 1] || self.element );
                        $panel.data('destroy.ui-tabs', true);
                    }
                    self.$panels = self.$panels.add( $panel );
                }
                // invalid tab href
                else
                    o.disabled.push(i + 1);
            });

            if (init) {

                // attach necessary classes for styling if not present
                $(this.element).hasClass(o.navClass) || $(this.element).addClass(o.navClass);
                this.$panels.each(function() {
                    var $this = $(this);
                    $this.hasClass(o.panelClass) || $this.addClass(o.panelClass);
                });

                // disabled tabs
                for (var i = 0, index; index = o.disabled[i]; i++)
                    this.disable(index);

                // Try to retrieve selected tab:
                // 1. from fragment identifier in url if present
                // 2. from cookie
                // 3. from selected class attribute on <li>
                // 4. otherwise use given "selected" option
                // 5. check if tab is disabled
                this.$tabs.each(function(i, a) {
                    if (location.hash) {
                        if (a.hash == location.hash) {
                            o.selected = i;
                            // prevent page scroll to fragment
                            //if (($.browser.msie || $.browser.opera) && !o.remote) {
                            if ($.browser.msie || $.browser.opera) {
                                var $toShow = $(location.hash), toShowId = $toShow.attr('id');
                                $toShow.attr('id', '');
                                setTimeout(function() {
                                    $toShow.attr('id', toShowId); // restore id
                                }, 500);
                            }
                            scrollTo(0, 0);
                            return false; // break
                        }
                    } else if (o.cookie) {
                        var index = parseInt($.cookie('ui-tabs' + $.data(self.element)),10);
                        if (index && self.$tabs[index]) {
                            o.selected = index;
                            return false; // break
                        }
                    } else if ( self.$lis.eq(i).hasClass(o.selectedClass) ) {
                        o.selected = i;
                        return false; // break
                    }
                });
                var n = this.$lis.length;
                while (this.$lis.eq(o.selected).hasClass(o.disabledClass) && n) {
                    o.selected = ++o.selected < this.$lis.length ? o.selected : 0;
                    n--;
                }
                if (!n) // all tabs disabled, set option unselect to true
                    o.unselect = true;

                // highlight selected tab
                this.$panels.addClass(o.hideClass);
                this.$lis.removeClass(o.selectedClass);
                if (!o.unselect) {
                    this.$panels.eq(o.selected).show().removeClass(o.hideClass); // use show and remove class to show in any case no matter how it has been hidden before
                    this.$lis.eq(o.selected).addClass(o.selectedClass);
                }

                // load if remote tab
                var href = !o.unselect && $.data(this.$tabs[o.selected], 'load.ui-tabs');
                if (href)
                    this.load(o.selected, href);

                // disable click if event is configured to something else
                if (!(/^click/).test(o.event))
                    this.$tabs.bind('click', function(e) { e.preventDefault(); });

            }

            var hideFx, showFx, baseFx = { 'min-width': 0, duration: 1 }, baseDuration = 'normal';
            if (o.fx && o.fx.constructor == Array)
                hideFx = o.fx[0] || baseFx, showFx = o.fx[1] || baseFx;
            else
                hideFx = showFx = o.fx || baseFx;

            // reset some styles to maintain print style sheets etc.
            var resetCSS = { display: '', overflow: '', height: '' };
            if (!$.browser.msie) // not in IE to prevent ClearType font issue
                resetCSS.opacity = '';

            // Hide a tab, animation prevents browser scrolling to fragment,
            // $show is optional.
            function hideTab(clicked, $hide, $show) {
                $hide.animate(hideFx, hideFx.duration || baseDuration, function() { //
                    $hide.addClass(o.hideClass).css(resetCSS); // maintain flexible height and accessibility in print etc.
                    if ($.browser.msie && hideFx.opacity)
                        $hide[0].style.filter = '';
                    if ($show)
                        showTab(clicked, $show, $hide);
                });
            }

            // Show a tab, animation prevents browser scrolling to fragment,
            // $hide is optional.
            function showTab(clicked, $show, $hide) {
                if (showFx === baseFx)
                    $show.css('display', 'block'); // prevent occasionally occuring flicker in Firefox cause by gap between showing and hiding the tab panels
                $show.animate(showFx, showFx.duration || baseDuration, function() {
                    $show.removeClass(o.hideClass).css(resetCSS); // maintain flexible height and accessibility in print etc.
                    if ($.browser.msie && showFx.opacity)
                        $show[0].style.filter = '';

                    // callback
                    $(self.element).triggerHandler("show.ui-tabs", [self.ui(clicked, $show[0])]);

                });
            }

            // switch a tab
            function switchTab(clicked, $li, $hide, $show) {
                /*if (o.bookmarkable && trueClick) { // add to history only if true click occured, not a triggered click
                    $.ajaxHistory.update(clicked.hash);
                }*/
                $li.addClass(o.selectedClass)
                    .siblings().removeClass(o.selectedClass);
                hideTab(clicked, $hide, $show);
            }

            // attach tab event handler, unbind to avoid duplicates from former tabifying...
            this.$tabs.unbind(o.event).bind(o.event, function() {

                //var trueClick = e.clientX; // add to history only if true click occured, not a triggered click
                var $li = $(this).parents('li:eq(0)'),
                    $hide = self.$panels.filter(':visible'),
                    $show = $(this.hash);

                // If tab is already selected and not unselectable or tab disabled or click callback returns false stop here.
                // Check if click handler returns false last so that it is not executed for a disabled tab!
                if (($li.hasClass(o.selectedClass) && !o.unselect) || $li.hasClass(o.disabledClass)
                    || $(self.element).triggerHandler("select.ui-tabs", [self.ui(this, $show[0])]) === false) {
                    this.blur();
                    return false;
                }

                self.options.selected = self.$tabs.index(this);

                // if tab may be closed
                if (o.unselect) {
                    if ($li.hasClass(o.selectedClass)) {
                        self.options.selected = null;
                        $li.removeClass(o.selectedClass);
                        self.$panels.stop();
                        hideTab(this, $hide);
                        this.blur();
                        return false;
                    } else if (!$hide.length) {
                        self.$panels.stop();
                        var a = this;
                        self.load(self.$tabs.index(this), function() {
                            $li.addClass(o.selectedClass).addClass(o.unselectClass);
                            showTab(a, $show);
                        });
                        this.blur();
                        return false;
                    }
                }

                if (o.cookie)
                    $.cookie('ui-tabs' + $.data(self.element), self.options.selected, o.cookie);

                // stop possibly running animations
                self.$panels.stop();

                // show new tab
                if ($show.length) {

                    // prevent scrollbar scrolling to 0 and than back in IE7, happens only if bookmarking/history is enabled
                    /*if ($.browser.msie && o.bookmarkable) {
                        var showId = this.hash.replace('#', '');
                        $show.attr('id', '');
                        setTimeout(function() {
                            $show.attr('id', showId); // restore id
                        }, 0);
                    }*/

                    var a = this;
                    self.load(self.$tabs.index(this), function() {
                        switchTab(a, $li, $hide, $show);
                    });

                    // Set scrollbar to saved position - need to use timeout with 0 to prevent browser scroll to target of hash
                    /*var scrollX = window.pageXOffset || document.documentElement && document.documentElement.scrollLeft || document.body.scrollLeft || 0;
                    var scrollY = window.pageYOffset || document.documentElement && document.documentElement.scrollTop || document.body.scrollTop || 0;
                    setTimeout(function() {
                        scrollTo(scrollX, scrollY);
                    }, 0);*/

                } else
                    throw 'jQuery UI Tabs: Mismatching fragment identifier.';

                // Prevent IE from keeping other link focussed when using the back button
                // and remove dotted border from clicked link. This is controlled in modern
                // browsers via CSS, also blur removes focus from address bar in Firefox
                // which can become a usability and annoying problem with tabsRotate.
                if ($.browser.msie)
                    this.blur();

                //return o.bookmarkable && !!trueClick; // convert trueClick == undefined to Boolean required in IE
                return false;

            });

        },
        add: function(url, label, index) {
            if (url && label) {
                index = index || this.$tabs.length; // append by default

                var o = this.options;
                var $li = $(o.tabTemplate.replace(/#\{href\}/, url).replace(/#\{label\}/, label));
                $li.data('destroy.ui-tabs', true);

                var id = url.indexOf('#') == 0 ? url.replace('#', '') : this.tabId( $('a:first-child', $li)[0] );

                // try to find an existing element before creating a new one
                var $panel = $('#' + id);
                if (!$panel.length) {
                    $panel = $(o.panelTemplate).attr('id', id)
                        .addClass(o.panelClass).addClass(o.hideClass);
                    $panel.data('destroy.ui-tabs', true);
                }
                if (index >= this.$lis.length) {
                    $li.appendTo(this.element);
                    $panel.appendTo(this.element.parentNode);
                } else {
                    $li.insertBefore(this.$lis[index]);
                    $panel.insertBefore(this.$panels[index]);
                }

                this.tabify();

                if (this.$tabs.length == 1) {
                     $li.addClass(o.selectedClass);
                     $panel.removeClass(o.hideClass);
                     var href = $.data(this.$tabs[0], 'load.ui-tabs');
                     if (href)
                         this.load(index, href);
                }

                // callback
                $(this.element).triggerHandler("add.ui-tabs",
                    [this.ui(this.$tabs[index], this.$panels[index])]
                );

            } else
                throw 'jQuery UI Tabs: Not enough arguments to add tab.';
        },
        remove: function(index) {
            if (index && index.constructor == Number) {
                var o = this.options, $li = this.$lis.eq(index).remove(),
                    $panel = this.$panels.eq(index).remove();

                // If selected tab was removed focus tab to the right or
                // tab to the left if last tab was removed.
                if ($li.hasClass(o.selectedClass) && this.$tabs.length > 1)
                    this.click(index + (index < this.$tabs.length ? 1 : -1));
                this.tabify();

                // callback
                $(this.element).triggerHandler("remove.ui-tabs",
                    [this.ui($li.find('a')[0], $panel[0])]
                );

            }
        },
        enable: function(index) {
            var self = this, o = this.options, $li = this.$lis.eq(index);
            $li.removeClass(o.disabledClass);
            if ($.browser.safari) { // fix disappearing tab (that used opacity indicating disabling) after enabling in Safari 2...
                $li.css('display', 'inline-block');
                setTimeout(function() {
                    $li.css('display', 'block');
                }, 0);
            }

            o.disabled = $.map(this.$lis.filter('.' + o.disabledClass),
                function(n, i) { return self.$lis.index(n); } );

            // callback
            $(this.element).triggerHandler("enable.ui-tabs",
                [this.ui(this.$tabs[index], this.$panels[index])]
            );

        },
        disable: function(index) {
            var self = this, o = this.options;
            this.$lis.eq(index).addClass(o.disabledClass);

            o.disabled = $.map(this.$lis.filter('.' + o.disabledClass),
                function(n, i) { return self.$lis.index(n); } );

            // callback
            $(this.element).triggerHandler("disable.ui-tabs",
                [this.ui(this.$tabs[index], this.$panels[index])]
            );

        },
        select: function(index) {
            if (typeof index == 'string')
                index = this.$tabs.index( this.$tabs.filter('[href$=' + index + ']')[0] );
            this.$tabs.eq(index).trigger(this.options.event);
        },
        load: function(index, callback) { // callback is for internal usage only
            var self = this, o = this.options,
                $a = this.$tabs.eq(index), a = $a[0];

            var url = $a.data('load.ui-tabs');

            // no remote - just finish with callback
            if (!url) {
                typeof callback == 'function' && callback();
                return;
            }

            // load remote from here on
            if (o.spinner) {
                var $span = $('span', a), label = $span.html();
                $span.html('<em>' + o.spinner + '</em>');
            }
            var finish = function() {
                self.$tabs.filter('.' + o.loadingClass).each(function() {
                    $(this).removeClass(o.loadingClass);
                    if (o.spinner)
                        $('span', this).html(label);
                });
                self.xhr = null;
            };
            var ajaxOptions = $.extend({}, o.ajaxOptions, {
                url: url,
                success: function(r, s) {
                    $(a.hash).html(r);
                    finish();
                    // This callback is required because the switch has to take
                    // place after loading has completed.
                    typeof callback == 'function' && callback();

                    if (o.cache)
                        $.removeData(a, 'load.ui-tabs'); // if loaded once do not load them again

                    // callback
                    $(self.element).triggerHandler("load.ui-tabs",
                        [self.ui(self.$tabs[index], self.$panels[index])]
                    );

                    o.ajaxOptions.success && o.ajaxOptions.success(r, s);
                }
            });
            if (this.xhr) {
                // terminate pending requests from other tabs and restore tab label
                this.xhr.abort();
                finish();
            }
            $a.addClass(o.loadingClass);
            setTimeout(function() { // timeout is again required in IE, "wait" for id being restored
                self.xhr = $.ajax(ajaxOptions);
            }, 0);

        },
        url: function(index, url) {
            this.$tabs.eq(index).data('load.ui-tabs', url);
        },
        destroy: function() {
            var o = this.options;
            $(this.element).unbind('.ui-tabs')
                .removeClass(o.navClass).removeData('ui-tabs');
            this.$tabs.each(function() {
                var href = $.data(this, 'href.ui-tabs');
                if (href)
                    this.href = href;
                $(this).unbind('.ui-tabs')
                    .removeData('href.ui-tabs').removeData('load.ui-tabs');
            });
            this.$lis.add(this.$panels).each(function() {
                if ($.data(this, 'destroy.ui-tabs'))
                    $(this).remove();
                else
                    $(this).removeClass([o.selectedClass, o.unselectClass,
                        o.disabledClass, o.panelClass, o.hideClass].join(' '));
            });
        }
    });

})(jQuery);

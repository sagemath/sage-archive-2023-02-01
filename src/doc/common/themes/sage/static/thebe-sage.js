$(function() {
    if (window.location.protocol.startsWith('http')) {
        $('<button id="thebe-activate">Activate</button>')
            .css({position: 'absolute', right: 0})
            .click(function() {
                new Thebe({
                     tmpnb_mode: false,
                     load_css: false,
                     url: window.location.origin,
                     kernel_name: "sagemath",
                     selector: "pre:contains('sage: ')"
                });
            })
            .prependTo('div.body');
    }
});

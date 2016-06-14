$(function() {

    function setupThebeSage(thebeOptions) {
        return new Thebe(thebeOptions);
    }

    if (window.location.protocol.startsWith('http')) {
        var selector = "pre:contains('sage: ')";
        if ($(selector).length > 0) {
            $('<button id="thebe-activate">Activate</button>')
                .css({position: 'absolute', right: 0})
                .click(function() {
                    setupThebeSage({
                        tmpnb_mode: false,
                        load_css: false,
                        url: window.location.origin,
                        kernel_name: "sagemath",
                        selector: selector
                    });
                    $(this).attr('disabled', 'disabled');
                })
                .prependTo('div.body');
        }
    }
});

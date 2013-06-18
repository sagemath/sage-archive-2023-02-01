function generic_callback(status, response_text) {
    /* do nothing */
}

function async_request(url, callback, postvars) {
    var settings = {url : url,
                    async : true,
                    cache : false,
                    dataType: "text"};

    if($.isFunction(callback)) {
        settings['error'] = function (XMLHttpRequest, textStatus, errorThrown) {
            callback("failure", errorThrown);
        };
        settings['success'] = function (data, textStatus) {
            callback("success", data);
        }
    }

    if(postvars != null) {
        settings['type'] = "POST";
        settings['data'] = postvars;
    } else {
        settings['type'] = "GET";
    }

    $.ajax(settings);
}

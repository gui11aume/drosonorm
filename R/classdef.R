setClass ("marrayPlatform",
    representation = representation(
        name = "character",
        release = "character",
        mapping = "data.frame",
        spikes = "character",
        random = "character"));

# Overload '$'.
getGeneric('$');
setMethod('$',
    signature(x = "marrayPlatform"),
    definition = function(x, name) {
        return (slot(x, name));
    });

(function() {
  var eventEmitter = new EventEmitter();

  window.sample_defaults = {
    addListener: function(event, listener) {
      eventEmitter.addListener(event, listener);
    },
    // width: 320,
    // height: 240,
    paused: false,
    wireframe: false,
    // current_normal_map: "normal_map_tile.jpg",
    // normal_maps: [
    //   "normal_map_face.PNG",
    //   "normal_map_circle.jpg",
    //   "normal_map_tile.jpg"
    // ]
  };

  window.samples = {};
  function createSample($el) {
    var index = $el.data("sample");
    var instance = window.samples[index].initialize($el[0]);
    $el.data("instance", instance);
    return instance;
  };

  function runCurrentSample(currentSlide) {
    $(currentSlide).find("[data-sample]").each(function() {
      var instance = createSample($(this));
      if(instance) instance.active = true;
    });
  };

  function activateCurrentSample(currentSlide) {
    $(currentSlide).find("[data-sample]").each(function() {
      console.log("Activating sample");
      var instance = $(this).data("instance");
      if(instance) instance.active = true;
    });
  }

  function initializeOnLoad() {
    // runCurrentSample($("section.present"));
    // Activate appropriate sample on slide change.
    Reveal.addEventListener('slidechanged', function(event) {
      // Clear all slides
      $("[data-sample]").each(function() {
        var instance = $(this).data("instance");
        if(instance) instance.active = false;
      });

      var currentSlide = event.currentSlide;
      // runCurrentSample(currentSlide);
      activateCurrentSample(currentSlide);
    });

    eventEmitter.emitEvent("initialized");
  }

  function initializeAll() {
    $("[data-sample]").each(function() {
      var instance = createSample($(this));
      if(instance) instance.active = false;
    });
    Reveal.addEventListener('slidechanged', function(event) {
      // Clear all slides
      $("[data-sample]").each(function() {
        var instance = $(this).data("instance");
        if(instance) instance.active = false;
      });

      var currentSlide = event.currentSlide;
      // runCurrentSample(currentSlide);
      activateCurrentSample(currentSlide);
    });

    eventEmitter.emitEvent("initialized");
  }

  head.js(
    // "js/samples/simulation.js",
    "js/samples/star_dipole.js",
    // "js/samples/twisted_dipole.js",
    // "js/samples/radiation.js",
    // "js/samples/pair_creation.js",
    // initializeOnLoad,
    initializeAll);
})();

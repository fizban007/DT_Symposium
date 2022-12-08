(function() {

    window.samples.twisted_dipole = {

	initialize: function(canvas) {
	    var scene = new THREE.Scene();
	    scene.background = new THREE.Color( 0x222222 );
	    // THREE.Object3D.DefaultUp.set(0.5,0.0,0.8);
	    var width = 1000;
	    var height = 720;
	    var camera = new THREE.PerspectiveCamera( 30, width / height, 1, 1000 );
	    camera.position.z = 0;
	    camera.position.y = 25;
	    camera.up.set(0,0,1);
	    camera.lookAt( new THREE.Vector3(0,0,0));

	    var radius = 1.0;
	    var geometry = new THREE.SphereGeometry( radius, 64, 64 );
	    var material = new THREE.MeshPhongMaterial({ color: 0xaaaaaa });
	    // var material = new THREE.MeshBasicMaterial({ color: 0xaaaaaa });
	    var mesh = new THREE.Mesh( geometry, material );
	    scene.add( mesh );

	    function add_dipole_line(theta_i, phi, smoothness, color) {
		var s2i = Math.sin(theta_i) * Math.sin(theta_i);

		var angleRange = Math.PI - 2 * theta_i; // get the angle between vectors
		var angleDelta = angleRange / (smoothness - 1); // increment

		var g = new THREE.Geometry();
		for (var i = 0; i < smoothness; i++) {
		    var th = i * angleDelta + theta_i;
		    var r = Math.sin(th) * Math.sin(th) * radius / s2i;
		    var p = new THREE.Vector3(r*Math.sin(th)*Math.cos(phi),
					      r*Math.sin(th)*Math.sin(phi),
					      r*Math.cos(th));
		    g.vertices.push(p)  // this is the key operation
		}

		var arc = new THREE.Line(g, new THREE.LineBasicMaterial({
		    color: color,
		    linewidth: 1.5
		}));
	    	return arc;
	    }

	    var field_lines = new THREE.Group();
	    var smoothness = 70;
	    for (var i = 0; i < 5; i++) {
		var theta = 0.25 + i * 0.12;
		for (var j = 0; j < 2; j++) {
		    var phi = 0.2 + 0.5 * Math.PI + j * 0.05;
		    var line = add_dipole_line(theta, phi, smoothness, "limegreen");
		    field_lines.add(line);
		}
	    }
	    scene.add(field_lines);

	    var twist = [];
	    var twist_dt = 0.009;
	    for (var i = 0; i < smoothness; i++) {
		twist.push((2.0 * i / (smoothness - 1) - 1.0) * twist_dt);
	    }

	    var directionalLight = new THREE.DirectionalLight ( 0xffffffff );
	    directionalLight.position.set( 0, 7, 0);
	    scene.add( directionalLight );

	    var renderer = new THREE.WebGLRenderer({canvas: canvas, antialias: true});
	    renderer.setSize( width , height );
	    var controls = new THREE.OrbitControls(camera, renderer.domElement);
	    controls.enableKeys = false;

	    var instance = { active: false, finished: false };

	    var numFrames = 100, frame = 0, paused = true, finished = false;
	    function animate() {
		requestAnimationFrame( animate, canvas );
		if(!instance.active || sample_defaults.paused) return;

		// field_lines.rotateOnWorldAxis(new THREE.Vector3(0,0,1), 0.003);
		field_lines.traverse( function(node) {
		    if (node instanceof THREE.Line
		        && !paused) {
			// console.log("found one");
			// console.log(node.geometry.vertices.length);
			for (var i = 0; i < node.geometry.vertices.length; i++) {
			    node.geometry.vertices[i].applyAxisAngle(new THREE.Vector3(0,0,1),
								     -twist[i]-0.002);
			}
			node.geometry.verticesNeedUpdate = true;
		    }
		});
		if (!paused) {
		    frame += 1;
		    if (frame > numFrames) {
			paused = true;
			Reveal.removeKeyBinding(39);
			instance.finished = true;
		    }
		}
		console.log(frame);

		renderer.render( scene, camera );
	    }

	    // renderer.domElement.onkeydown = function(event) {
	    instance.onkeydown = function(event) {
		var keyCode = event.code || event.keyCode || 0;
		// Respond to either q or down arrow
		if (keyCode === 'KeyQ' ||
		    event.keyCode === 39) {
		    if (!instance.finished)
			paused = !paused;
		}
	    };
	    // renderer.domElement.addEventListener("mousedown", onClick, false);
	    // document.addEventListener("mousedown", onClick, false);
	    // function onClick(event) {
	    // 	console.log("clicked");
	    // 	paused = !paused;
	    // };

	    animate();
	    return instance;
	}
    };
})();

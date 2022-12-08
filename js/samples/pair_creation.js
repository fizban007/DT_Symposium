(function() {

    window.samples.pair_creation = {

	initialize: function(canvas) {
	    var scene = new THREE.Scene();
	    scene.background = new THREE.Color( 0x222222 );
	    // THREE.Object3D.DefaultUp.set(0.5,0.0,0.8);
	    var width = 1000;
	    var height = 600;
	    // var camera = new THREE.PerspectiveCamera( 30, width / height, 1, 1000 );
	    var camera = new THREE.OrthographicCamera( width / - 100, width / 100,
						       height / 100, height / - 100, 1, 1000 );
	    camera.position.z = 0;
	    camera.position.y = 3;
	    camera.position.x = -7;
	    camera.up.set(0,0,1);
	    camera.lookAt( new THREE.Vector3(-7,0,0));

	    var photon_arrow = new THREE.TextureLoader().load('figures/photon-arrow.png');
	    var fieldline_txt = new THREE.TextureLoader().load('figures/fieldline.png');

	    function add_field_line(x, length, smoothness, color) {
		var g = new THREE.Geometry();
		for (var i = 0; i < smoothness; i++) {
		    var p = new THREE.Vector3(x, 0, i * length/(smoothness - 1));
		    g.vertices.push(p)  // this is the key operation
		}

		var fline = new MeshLine();
		fline.setGeometry(g, function(p) { return 0.004; });
		var arc = new THREE.Mesh(fline.geometry, new MeshLineMaterial({
		    // color: color,
		    useMap: true,
		    map: fieldline_txt,
		    transparent: true,
		}));
	    	return arc;
	    }

	    function add_animated_photon_path(start, end, smoothness, color) {
		var g = new THREE.Geometry();
		for (var i = 0; i < smoothness; i++) {
		    var r = r_i + (r_f - r_i) * i / (smoothness - 1);
		    var p = new THREE.Vector3(-r, 0, 0.08*Math.sin(2*Math.PI*r/(0.15*(r_f - r_i))));
		    g.vertices.push(p);
		}
		var phline = new MeshLine();
		phline.setGeometry(g, function(p) { return 0.005; });
		var ph = new THREE.Mesh(phline.geometry, new MeshLineMaterial({
		    useMap: true,
		    map: photon_arrow,
		    transparent: true,
		    opacity: 0.9,
		    // color: color,
		    dashArray: 1.5,     // always has to be the double of the line
		    dashOffset: 0,    // start the dash at zero
		    dashRatio: 0.7,  // visible length range min: 0.99, max: 0.5
		}));
		ph.rotateY(theta);
		ph.rotateZ(phi);
		return ph;
	    }

	    function add_photon_path(start, end, smoothness, color) {
		var g = new THREE.Geometry();
		for (var i = 0; i < smoothness; i++) {
		    var x = i / (smoothness - 1);
		    // var r = r_i + (r_f - r_i) * i / (smoothness - 1);
		    var p = start + x * (end - start);
		    g.vertices.push(p);
		}
		var ph = new THREE.Line(g, new THREE.LineBasicMaterial({
		    color: color
		}));
		var origin = g.vertices[smoothness - 1];
		var length = 0.2*(r_f - r_i);
		var dir = new THREE.Vector3(-1, 0, 0);
		var arr = new THREE.ArrowHelper(dir, origin, length, color.getHex(), 0.99*length, 0.5*length);
		ph.add(arr);
		ph.rotateY(theta);
		ph.rotateZ(phi);
		return ph;
	    }

	    var photons = new THREE.Group();
	    var ph_smooth = 100;
	    // var ph = add_photon_path(0, 0, 0.5, 3.5, ph_smooth, new THREE.Color("red"));
	    // photons.add(ph);
	    for (var i = 0; i < 6; i++) {
	    	var ph = add_animated_photon_path(-1.25 + 0.5*i, 0.0, 1.0, 3.5, ph_smooth, new THREE.Color("red"));
	    	photons.add(ph);
	    }
	    photons.visible = false;
	    scene.add(photons);

	    var electron_geom = new THREE.SphereGeometry(0.1, 32, 32);
	    electron_geom.translate(-2.8, 0, 2.75);
	    var electron = new THREE.Mesh(electron_geom, material);
	    electron.visible = false;
	    scene.add(electron);

	    var gamma_ray = add_photon_path(0.2, 0.0, 4.0, 5.7, ph_smooth, new THREE.Color("lightblue"));
	    gamma_ray.translateX(0.2);
	    gamma_ray.translateZ(2.2);
	    gamma_ray.visible = false;
	    scene.add(gamma_ray);

	    var field_lines = new THREE.Group();
	    var smoothness = 70;
	    for (var i = 0; i < 5; i++) {
		var theta = 0.25 + i * 0.12;
		for (var j = 0; j < 2; j++) {
		    var phi = 0.2 + 0.9 * Math.PI + j * 0.05;
		    var line = add_dipole_line(theta, phi, smoothness, new THREE.Color("limegreen"));
		    field_lines.add(line);
		}
	    }
	    scene.add(field_lines);

	    var directionalLight = new THREE.DirectionalLight ( 0xffffffff );
	    directionalLight.position.set( 0, 7, 0);
	    scene.add( directionalLight );

	    var renderer = new THREE.WebGLRenderer({canvas: canvas, antialias: true});
	    renderer.setSize( width , height );
	    // var controls = new THREE.OrbitControls(camera, renderer.domElement);
	    // controls.enableKeys = false;

	    var instance = { active: false };

	    var stage = 0;

	    function animate() {
		requestAnimationFrame( animate, canvas );
		if(!instance.active || sample_defaults.paused) return;

		if (stage == 1) {
		    photons.visible = true;
		} else if (stage == 2) {
		    electron.visible = true;
		} else if (stage == 3) {
		    gamma_ray.visible = true;
		}

		if (photons.visible == true && gamma_ray.visible == false) {
		    photons.traverse(function(node){
			if (node instanceof THREE.Mesh
			    && node.material instanceof MeshLineMaterial) {
			    // Check if the dash is out to stop animate it.
			    if (node.material.uniforms.dashOffset.value >= -1)
				node.material.uniforms.dashOffset.value -= 0.01;
				// node.material.uniforms.dashOffset.value = 0;
			    // console.log(node.material.uniforms.dashOffset.value);

			    // Decrement the dashOffset value to animate the path with the dash.
			    node.material.needsUpdate = true;
			}
		    });
		}

		renderer.render( scene, camera );
	    }

	    instance.onkeydown = function(event) {
		var keyCode = event.code || event.keyCode || 0;
		// Respond to either q or down arrow
		if (keyCode === 'KeyQ' ||
		    event.keyCode === 39) {
		    stage += 1;
		    if (stage >= 4) {
			// photons.visible = false;
			// electron.visible = false;
			// gamma_ray.visible = false;
			// stage = 0;
			Reveal.next();
		    }
		}
	    };

	    animate();

	    return instance;
	}
    };
})();

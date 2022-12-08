(function() {

    window.samples.simulation = {

	initialize: function(canvas) {
	    var folder = "data/magnetar_talk/";
	    var folder_lines = "data/magnetar_talk/lines/"
	    var normalization = 1;
	    var shift = 0;
	    var total_steps = 140000;
	    var paused = true;
	    var rot_paused = true;

	    var scene = new THREE.Scene();
	    scene.background = new THREE.Color( 0x222222 );
	    // THREE.Object3D.DefaultUp.set(0.5,0.0,0.8);
	    var width = window.innerWidth;
	    var height = window.innerHeight;
	    const aspect = width / height;
	    // var camera = new THREE.PerspectiveCamera( 30, width / height, 1, 1000 );
	    var camera = new THREE.PerspectiveCamera( 75, aspect, 0.1, 1000 );
	    var cam_radius = 10;
	    camera.position.y = -cam_radius;
	    // camera.position.x = 7;
	    camera.lookAt([0, 0, 0]);
	    camera.up = new THREE.Vector3(0, 0, 1);
	    camera.updateProjectionMatrix();
	    var controls = new THREE.OrbitControls( camera );

	    controls.update();
	    controls.enableKeys = false;

	    var renderer = new THREE.WebGLRenderer({
		canvas: canvas,
		antialias: true,
		preserveDrawingBuffer: true,
	    });
	    renderer.setSize( width , height );

	    // ------------------------------------------------
	    // FUN STARTS HERE
	    // ------------------------------------------------
	    var particleSystem, uniforms, pGeometry;
	    var particles = 1000000;
	    var vert_shader, frag_shader;
	    var pos_e = new Float32Array(particles * 3), pos_p = new Float32Array(particles * 3);
	    var geometry_e = new THREE.BufferGeometry();
	    var geometry_p = new THREE.BufferGeometry();
	    geometry_e.addAttribute( 'position', new THREE.BufferAttribute(pos_e, 3) );
	    geometry_p.addAttribute( 'position', new THREE.BufferAttribute(pos_p, 3) );

	    // Create a Sphere Mesh with basic material
	    var geometry = new THREE.SphereGeometry( 1.03, 32, 32 );
	    // var material = new THREE.MeshBasicMaterial( { color: "#888888",
	    //                                               depthTest: true,
	    //                                               depthWrite: true,
	    //                                               transparent: true } );
	    var material = new THREE.MeshPhongMaterial({ color: 0x888888 });
	    var sphere = new THREE.Mesh( geometry, material );
	    sphere.matrixAutoUpdate = false;
	    // Add cube to Scene
	    scene.add( sphere );

	    var instance = { active: false };

	    // Load the vertex and fragment shaders from external file
	    ShaderLoader("shaders/particle.vert", "shaders/particle.frag", function (vs, fs) {
		start(vs, fs);
	    });

	    function start(vs, fs) {
		var uniforms_e = {
		    texture:  { value: new THREE.TextureLoader().load( "textures/sprites/spark1.png" ) },
		    uCol: new THREE.Uniform(new THREE.Color(0.02, 0.02, 1.0))
		};
		var uniforms_p = {
		    texture:  { value: new THREE.TextureLoader().load( "textures/sprites/spark1.png" ) },
		    uCol: new THREE.Uniform(new THREE.Color(1.0, 0.02, 0.02))
		};
		var shaderMaterial_e = new THREE.ShaderMaterial( {
		    uniforms:       uniforms_e,
		    vertexShader:   vs,
		    fragmentShader: fs,
		    blending:       THREE.CustomBlending,
		    blendEquation:  THREE.AddEquation,
		    blendSrc:       THREE.SrcAlphaFactor,
		    blendDst:       THREE.OneFactor,
		    depthWrite:     false,
		    transparent:    true,
		    vertexColors:   true
		});
		var shaderMaterial_p = new THREE.ShaderMaterial( {
		    uniforms:       uniforms_p,
		    vertexShader:   vs,
		    fragmentShader: fs,
		    blending:       THREE.CustomBlending,
		    blendEquation:  THREE.AddEquation,
		    blendSrc:       THREE.SrcAlphaFactor,
		    blendDst:       THREE.OneFactor,
		    depthWrite:     false,
		    transparent:    true,
		    vertexColors:   true
		});
		var radius = 5;

		var electrons = new THREE.Points(geometry_e, shaderMaterial_e);
		var positrons = new THREE.Points(geometry_p, shaderMaterial_p);

		scene.add(electrons);
		scene.add(positrons);

		// Read the files into array buffers
		var loader = new THREE.FileLoader();
		loader.setResponseType('arraybuffer');
		var fieldline_txt = new THREE.TextureLoader().load('figures/fieldline.png');
		var line_mat = new MeshLineMaterial({
		    useMap: true,
		    map: fieldline_txt,
		    // color: new THREE.Color(0x6aed5a),
		    opacity: 0.9,
		    depthWrite: false,
		    // depthTest: false,
		    transparent: true
		});
		var num_lines = 6;

		var prev_len_e = pos_e.length;
		var prev_len_p = pos_p.length;
		// Read an arraybuffer of particle positions
		function update_particles(data, species) {
		    // var floatView = new Float32Array(data);
		    var floatView = data;
		    var len = data.length;
		    if (species == 'e') {
			for (var i = 0; i < len/3; i++) {
			    pos_e[i * 3] = (floatView[i * 3] + shift) / normalization;
			    pos_e[i * 3 + 1] = (floatView[i * 3 + 1] + shift) / normalization;
			    pos_e[i * 3 + 2] = (floatView[i * 3 + 2] + shift) / normalization;
			}
			for (var i = len; i < prev_len_e; i++) {
			    pos_e[i] = 0.0;
			}
			prev_len_e = len;
			electrons.geometry.attributes.position.needsUpdate = true;
		    } else if (species == 'p') {
			for (var i = 0; i < len/3; i++) {
			    pos_p[i * 3] = (floatView[i * 3] + shift) / normalization;
			    pos_p[i * 3 + 1] = (floatView[i * 3 + 1] + shift) / normalization;
			    pos_p[i * 3 + 2] = (floatView[i * 3 + 2] + shift) / normalization;
			}
			for (var i = len; i < prev_len_p; i++) {
			    pos_p[i] = 0.0;
			}
			prev_len_p = len;
			positrons.geometry.attributes.position.needsUpdate = true;
		    }
		};

		function make_field_line(data, lines, line_meshes) {
		    var floatView = new Float32Array(data);
		    var geom = new THREE.Geometry();
		    var len = floatView.length;
		    for (var j = 0; j < len/3; j++) {
			var v = new THREE.Vector3(floatView[j*3], floatView[j*3+1], floatView[j*3+2]);
			geom.vertices.push(v);
		    }
		    var line = new MeshLine();
		    line.setGeometry(geom, function(p) {
			// return p * (1.0 - p);
			// return 0.1 * p * (1.0 - p);
			return 0.030;
		    });
		    var mesh = new THREE.Mesh(line.geometry, line_mat);
		    scene.add(mesh);
		    lines.push(line);
		    line_meshes.push(mesh);
		};

		function update_field_line(buffers, step, lines, line_meshes) {
		    for (var j = 0; j < num_lines; j += 1) {
			scene.remove(line_meshes[j]);
		    }
		    // line_meshes = [];
		    line_meshes.length = 0;
		    for (var j = 0; j < num_lines; j += 1) {
			var geom = new THREE.Geometry();
			// var g = lines[j].g;
			console.log("updating line", j, "at step", step, "length is", buffers[j][step].length);
			for (var k = 0; k < buffers[j][step].length/3; k += 1) {
			    var v = new THREE.Vector3(buffers[j][step][k*3],
						      buffers[j][step][k*3+1],
						      buffers[j][step][k*3+2]);
			    geom.vertices.push(v);
			}
			var line = new MeshLine();
			line.setGeometry(geom, function(p) {
			    return 0.030;
			});
			var newline = new THREE.Mesh(line.geometry, line_mat);
			scene.add(newline);
			line_meshes.push(newline);
			// line_meshes[j].geometry.attributes.position.needsUpdate = true;
			// line_meshes[j].geometry.elementsNeedUpdate = true;
		    }
		}

		var frame = 0;
		var startTime = Date.now();
		var pausedTime = Date.now();
		var rotation = 0.0;
		var lines = [];
		var line_meshes = [];
		loader.load(folder_lines + "line_0_000000", function(data) {make_field_line(data, lines, line_meshes);});
		loader.load(folder_lines + "line_1_000000", function(data) {make_field_line(data, lines, line_meshes);});
		loader.load(folder_lines + "line_2_000000", function(data) {make_field_line(data, lines, line_meshes);});
		loader.load(folder_lines + "line_3_000000", function(data) {make_field_line(data, lines, line_meshes);});
		loader.load(folder_lines + "line_4_000000", function(data) {make_field_line(data, lines, line_meshes);});
		loader.load(folder_lines + "line_5_000000", function(data) {make_field_line(data, lines, line_meshes);});
		// loader.load(folder_lines + "line_6_000000", function(data) {make_field_line(data, lines, line_meshes);});
		// loader.load(folder_lines + "line_7_000000", function(data) {make_field_line(data, lines, line_meshes);});
		// loader.load(folder_lines + "line_8_000000", function(data) {make_field_line(data, lines, line_meshes);});
		// loader.load(folder_lines + "line_9_000000", function(data) {make_field_line(data, lines, line_meshes);});

		var buffers_e = [];
		var buffers_p = [];
		// var buffer_lines = [[], [], [], [], [], [], [], [], [], []];
		var buffer_lines = [[], [], [], [], [], []];
		// var buffer_lines = [];

		var load_e_step = 500;
		var load_p_step = 500;
		var load_line_step = 500;
		var load_line_num = 0;
		function load_next_e() {
		    if (load_e_step >= total_steps) return;

		    console.log("loading electrons step ", load_e_step);
		    loader.load(folder + "pos_e_" + ("000" + load_e_step).slice(-6),
				function(data) {
				    var floatView = new Float32Array(data);
				    buffers_e.push(floatView);
				    load_e_step += 500;
				    load_next_e();
				});
		}
		function load_next_p() {
		    if (load_p_step >= total_steps) return;

		    console.log("loading positrons step ", load_p_step);
		    loader.load(folder + "pos_p_" + ("000" + load_p_step).slice(-6),
				function(data) {
				    var floatView = new Float32Array(data);
				    buffers_p.push(floatView);
				    load_p_step += 500;
				    load_next_p();
				});
		}
		function load_next_line() {
		    if (load_line_step >= total_steps) return;
		    console.log("loading line", load_line_num, "at step ", load_line_step);
		    loader.load(folder_lines + "line_" + load_line_num + "_" + ("000" + load_line_step).slice(-6),
				function(data) {
				    var floatView = new Float32Array(data);
				    buffer_lines[load_line_num].push(floatView);
				    if (load_line_num < num_lines - 1) {
					load_line_num += 1;
				    } else {
					load_line_num = 0;
					load_line_step += 500;
				    }
				    load_next_line();
				});
		}
		load_next_e();
		load_next_p();
		load_next_line();
		var directionalLight = new THREE.DirectionalLight ( 0xffffffff );
		directionalLight.position.set( 8, -7, 0);
		scene.add( directionalLight );
		var ambientLight = new THREE.AmbientLight( 0x404040 ); // soft white light
		scene.add( ambientLight );

		function onWindowResize() {
		    camera.aspect = window.innerWidth / window.innerHeight;
		    camera.updateProjectionMatrix();
		    renderer.setSize( window.innerWidth, window.innerHeight );
		}
		var then = 0;

		function animate(now) {
		    requestAnimationFrame( animate );
		    if(!instance.active || sample_defaults.paused) return;
		    now *= 0.001;
		    var delta = now - then;
		    // var time = Date.now();
		    // if (!paused && time - startTime > 1000/20) {
		    if (!paused && delta > 1.0/24.0) {
			// startTime = Date.now();
			if (frame < total_steps) {
			    frame += 500;
			    update_particles(buffers_e[frame/500], 'e');
			    update_particles(buffers_p[frame/500], 'p');
			    update_field_line(buffer_lines, frame/500, lines, line_meshes);
			    var time_text = document.getElementById("sim-time");
			    time_text.innerText = "t = " + (frame * 0.5 / 500);
			} else {
			    frame = 0;
			    for (var i = 0; i < pos_e.length; i++) {
				pos_e[i] = 0.0;
				pos_p[i] = 0.0;
			    }
			    paused = true;
			}
			console.log(frame);
			then = now;
		    }

		    // required if controls.enableDamping or controls.autoRotate are set to true
		    controls.update();
		    camera.updateProjectionMatrix();
		    render();
		    // stats.update();
		}

		function render() {
		    var time = (Date.now() - startTime) * 0.004;
		    startTime = Date.now();
		    // frame += 200;
		    positrons.visible = true;
		    electrons.visible = true;
		    // lines.forEach(function(l, i){
		    // 	l.visible = menu.fieldLines;
		    // });
		    // sphere.visible = menu.star;

		    if (!rot_paused) {
			cam_radius = Math.sqrt(camera.position.x * camera.position.x +
					       camera.position.y * camera.position.y);
			rotation = Math.atan(-camera.position.x/camera.position.y);
			camera.updateProjectionMatrix();
		    }
		    renderer.render( scene, camera );
		}

		document.body.onkeyup = function(event) {
		    var key = event.which || event.keyCode || 0;
		    if (key === 32) {
			paused = !paused;
			rot_paused = paused;
			if (!paused) {
			    startTime = Date.now();
			}
		    } else if (key === 82) {
			rot_paused = !rot_paused;
		    }
		};

		animate();
	    }

	    // This is a basic asyncronous shader loader for THREE.js.
	    function ShaderLoader(vertex_url, fragment_url, onLoad, onProgress, onError) {
		var vertex_loader = new THREE.FileLoader(THREE.DefaultLoadingManager);
		vertex_loader.setResponseType('text');
		vertex_loader.load(vertex_url, function (vertex_text) {
		    var fragment_loader = new THREE.FileLoader(THREE.DefaultLoadingManager);
		    fragment_loader.setResponseType('text');
		    fragment_loader.load(fragment_url, function (fragment_text) {
			onLoad(vertex_text, fragment_text);
		    });
		}, onProgress, onError);
	    }

	    return instance;
	}
    };
})();

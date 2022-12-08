(function() {

  window.samples.star_dipole = {

    initialize: function(canvas) {
      var scene = new THREE.Scene();
      scene.background = new THREE.Color( 0x222222 );
      // THREE.Object3D.DefaultUp.set(0.5,0.0,0.8);
      var width = window.innerWidth;
      var height = window.innerHeight;
      // var width = canvas.width;
      // var height = canvas.height;
      // const width = 1920;
      // const height = 1080;
      console.log(width, height);
      var camera = new THREE.PerspectiveCamera( 30, width / height, 1, 1000 );
      camera.position.z = 0;
      camera.position.y = 100;
      camera.up.set(0,0,1);
      camera.lookAt( new THREE.Vector3(0,0,0));

      var radius = 1.0;
      var geometry = new THREE.SphereGeometry( radius, 64, 64 );
      var material = new THREE.MeshPhongMaterial({ color: 0xaaaaaa });
      // var material = new THREE.MeshBasicMaterial({ color: 0xaaaaaa });
      var mesh = new THREE.Mesh( geometry, material );
      scene.add( mesh );

      var fieldline_txt = new THREE.TextureLoader().load('figures/fieldline.png');

      function add_dipole_line(theta_i, phi, smoothness, color) {
        var s2i = Math.sin(theta_i) * Math.sin(theta_i);

        var angleRange = Math.PI - 2 * theta_i; // get the angle between vectors
        var angleDelta = angleRange / (smoothness - 1); // increment

        var g = new THREE.Geometry();
        for (var i = 0; i < smoothness; i++) {
          var th = i * angleDelta + theta_i;
          var r = Math.sin(th) * Math.sin(th) * radius / s2i;
          // console.log(r, th);
          var p = new THREE.Vector3(r*Math.sin(th)*Math.cos(phi),
                                    r*Math.sin(th)*Math.sin(phi),
                                    r*Math.cos(th));
          g.vertices.push(p)  // this is the key operation
        }

        var fline = new MeshLine();
        fline.setGeometry(g, function(p) { return 0.3; });
        var arc = new THREE.Mesh(fline.geometry, new MeshLineMaterial({
          // color: color,
          useMap: true,
          map: fieldline_txt,
          transparent: true,
        }));
        return arc;
      }

      var field_lines = new THREE.Group();
      for (var i = 0; i < 3; i++) {
        var theta = 0.15 + i * 0.12;
        var nj = 8;
        for (var j = 0; j < nj; j++) {
          var phi = 0.3 + j * Math.PI * 2.0 / nj;
          var line = add_dipole_line(theta, phi, 80, new THREE.Color("limegreen"));
          field_lines.add(line);
        }
      }
      var g_cone = new THREE.ConeGeometry(3, 90, 32);
      var g_cone2 = new THREE.ConeGeometry(3, 90, 32);
      var mat_cone = new THREE.MeshPhongMaterial({
        color: 0xaaaaff,
        transparent: true,
        opacity: 0.6,
        emissive: 0xdddddd,
      });
      var cone = new THREE.Mesh(g_cone, mat_cone);
      var cone2 = new THREE.Mesh(g_cone2, mat_cone);
      cone.translateZ(-45.0);
      cone2.translateZ(45.0);
      cone.rotateX(0.5 * Math.PI);
      cone2.rotateX(-0.5 * Math.PI);
      // cone2.rotateY(Math.PI);
      field_lines.add(cone);
      field_lines.add(cone2);
      field_lines.rotateY(-0.8);
      scene.add(field_lines);

      var directionalLight = new THREE.DirectionalLight ( 0xffffffff );
      directionalLight.position.set( 0, 7, 0);
      scene.add( directionalLight );

      var light = new THREE.AmbientLight(0x404040); // soft white light
      scene.add(light);

      var renderer = new THREE.WebGLRenderer({canvas: canvas, antialias: true});
      renderer.setSize( width , height );
      var controls = new THREE.OrbitControls(camera, renderer.domElement);
      controls.enableKeys = false;

      var instance = { active: false };

      // console.log("Initialized star-dipole!");

      // canvas.addEventListener( 'resize', onCanvasResize, false );

      function resizeRendererToDisplaySize(renderer) {
        const canvas = renderer.domElement;
        const width = canvas.clientWidth;
        const height = canvas.clientHeight;
        const needResize = canvas.width !== width || canvas.height !== height;
        if (needResize) {
          renderer.setSize(width, height, false);
        }
        return needResize;
      }

      function onCanvasResize(){

        camera.aspect = canvas.width / canvas.height;
        camera.updateProjectionMatrix();

        renderer.setSize( canvas.width, canvas.height, false );

      }

      function animate() {
        if (resizeRendererToDisplaySize(renderer)) {
          const canvas = renderer.domElement;
          camera.aspect = canvas.clientWidth / canvas.clientHeight;
          camera.updateProjectionMatrix();
        }

        requestAnimationFrame( animate, canvas );
        if(!instance.active || sample_defaults.paused) return;

        field_lines.rotateOnWorldAxis(new THREE.Vector3(0,0,1), 0.01);

        renderer.render( scene, camera );
      }

      animate();
      return instance;
    }
  };
})();

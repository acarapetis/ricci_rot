'use strict';
var scene = new THREE.Scene();
var camera = new THREE.PerspectiveCamera( 75, window.innerWidth / window.innerHeight, 0.1, 1000 );
var controls = new THREE.TrackballControls( camera );

var renderer = new THREE.WebGLRenderer({ 
  alpha:      true,
  antialias:  true,
  preserveDrawingBuffer: true,
});
renderer.setSize( window.innerWidth, window.innerHeight );
document.body.appendChild( renderer.domElement );

controls.rotateSpeed = 1.0;
controls.zoomSpeed = 1.2;
controls.panSpeed = 0.8;

controls.noZoom = false;
controls.noPan = false;

controls.staticMoving = true;
controls.dynamicDampingFactor = 0.3;

controls.keys = [ 65, 83, 68 ];

controls.addEventListener( 'change', render );

var ambientLight = new THREE.AmbientLight( 0x000000 );
scene.add( ambientLight );

var lights = [];
var light_markers = [];
for (var i = 0; i < 1; i++) {
    lights[i] = new THREE.PointLight(0xffffff, 1, 0);
    lights[i].position.set( 0, 300, 0 );
    scene.add(lights[i]);
    light_markers[i] = new THREE.Mesh(
        new THREE.SphereGeometry(4,5,5),
        new THREE.MeshBasicMaterial({ color: 0xffff00 })
    );
    scene.add(light_markers[i]);
    lights[i].marker = light_markers[i];
}

var N = 801;
var W = 32;

var geometry = new RevolutionBufferGeometry(1, 64, 64, s => [
    //20 * Math.cos(Math.PI * s) + 2 * Math.cos(Math.PI * 10 * s), 
    //20 * Math.sin(Math.PI * s)
    20 * Math.cos(Math.PI*s),
    20 * Math.sin(Math.PI*s) + 2 * Math.sin(Math.PI*s*10)
]);

geometry.dynamic = true;

var material = new THREE.MeshPhongMaterial({ color: 0x2194CE });
var mesh = new THREE.Mesh( geometry, material );

scene.add( mesh );
camera.position.x = 0;
camera.position.y = -30;
camera.position.z = 70;
var frame = 0;

function render() {
	renderer.render( scene, camera );
}

function animate() {
    frame++;
    requestAnimationFrame( animate );
    controls.update();
    for (var light of lights) {
        var delta = new THREE.Vector3(
            20 * Math.random() - 10,
            20 * Math.random() - 10,
            20 * Math.random() - 10
        );
        light.position.add(delta)
                      .normalize()
                      .multiplyScalar(300);

        light.marker.position.set(light.position.x,
            light.position.y,
            light.position.z);
        light.marker.position.multiplyScalar(0.3);
    }
    render();
}
animate();

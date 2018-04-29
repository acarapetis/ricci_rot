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
lights[ 0 ] = new THREE.PointLight( 0xffffff, 1, 0 );
lights[ 1 ] = new THREE.PointLight( 0xffffff, 1, 0 );
lights[ 2 ] = new THREE.PointLight( 0xffffff, 1, 0 );

lights[ 0 ].position.set( 0, 200, 0 );
lights[ 1 ].position.set( 100, 200, 100 );
lights[ 2 ].position.set( - 100, - 200, - 100 );

scene.add( lights[ 0 ] );
scene.add( lights[ 1 ] );
scene.add( lights[ 2 ] );

var N = 801;
var W = 32;

var geometry = new RevolutionBufferGeometry(1, 32, 32, s => [
    20 * Math.cos(Math.PI * s) + 2 * Math.cos(Math.PI * 10 * s), 
    20 * Math.sin(Math.PI * s) - 2 * Math.sin(Math.PI * 10 * s)
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

  //torus.rotation.x += 0.01;
  //torus.rotation.y += 0.1;
  //render();
  //if (frame == 5) window.open(renderer.domElement.toDataURL("image/png"), 'DNA_Screen');
}
animate();
render();

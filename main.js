'use strict';
'import test.js';
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

var ambientLight = new THREE.AmbientLight( 0x333333 );
scene.add( ambientLight );

var lpos = [
    [100, 250, 100],
    [250, -100, 100],
    [100, 100, -250]
];
var lights = [];
for (var i = 0; i < lpos.length; i++) {
    lights[i] = new THREE.PointLight(0x666666, 1, 0);
    lights[i].position.set(lpos[i][0], lpos[i][1], lpos[i][2]);
    scene.add(lights[i]);
}

var N = Module._resolution() - 1;
var W = 32;

var geometry = new RevolutionBufferGeometry(1, W, N, s => [
    //20 * Math.cos(Math.PI * s) + 2 * Math.cos(Math.PI * 10 * s), 
    //20 * Math.sin(Math.PI * s)
    20 * Math.cos(Math.PI*s/N),
    20 * Math.sin(Math.PI*s/N) + 2 * Math.sin(Math.PI*s*10/N)
]);

var sphere = function(s) {
    return [ 20 * Math.cos(Math.PI*s), 20 * Math.sin(Math.PI*s) ];
};

geometry.dynamic = true;

var material = new THREE.MeshPhongMaterial({ color: 0xDDEEFF });
var mesh = new THREE.Mesh( geometry, material );

scene.add( mesh );
camera.position.x = 0;
camera.position.y = -30;
camera.position.z = 70;
var frame = 0;

var x_ary = new Float64Array(Module.HEAPF64.buffer,
                             Module._coord_x(),
                             N * 8);
var y_ary = new Float64Array(Module.HEAPF64.buffer,
                             Module._coord_y(),
                             N * 8);

function update_mesh() {
    geometry.updateRevFn(i => [
        20*x_ary[i], 20*y_ary[i]
    ]);
}

function reset_shape(c3, c5) {
    var result = Module._reset_shape(c3, c5);
    update_mesh();
}

function render() {
	renderer.render(scene, camera);
}

var old_t = false;
var ok = false;
function animate(t) {
    if (!old_t) old_t = t - 16;
    var dt = (t - old_t) * 0.000005;
    if (dt > 0.0001) dt = 0.0001;
    frame++;
    requestAnimationFrame( animate );
    controls.update();
    if (ok) Module._step(dt);
    if (!Module._make_xy()) {
        Module._undo_step();
        ok = false;
    }
    update_mesh();
    render();
    old_t = t;
}

function em_ready() {
    reset_shape(0.2,0.3);
    requestAnimationFrame(animate);
    ok = true;
}

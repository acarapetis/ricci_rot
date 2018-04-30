# ricci_rot
WebGL port of Robert Sinclair's ricci_rot

In [their 2005 paper](https://projecteuclid.org/euclid.em/1128371754), J. Hyam Rubenstein and Robert Sinclair proposed the visualization of Ricci flow via isometric embeddings of axially symmetric metrics, which have the nice property of remaining embeddable when acted upon by the flow. The paper was accompanied by Sinclair's ricci_rot program, which neatly demonstrated this visualization using OpenGL.

This is a (partial) adaptation of the ricci_rot demonstration to run in the browser, made possible by WebGL and Emscripten. The numerical simulation is exactly as in the original, while the visual and interactive components have been rewritten in Javascript using THREE.js.

## Building
Once you have installed a recent version of emscripten and configured your environment accordingly (e.g. via `source /path/to/emsdk_env.sh`), run `./build.sh` to compile `ricci_rot.c` to `ricci_rot.js`.

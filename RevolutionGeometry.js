function RevolutionBufferGeometry(height, radialSegments, heightSegments, revFn) { 

	THREE.BufferGeometry.call( this );

	this.type = 'RevolutionBufferGeometry';

	this.parameters = {
		height: height,
		radialSegments: radialSegments,
		heightSegments: heightSegments,
        revFn: revFn,
	};

	var scope = this;

	height = height || 1;

	this.radialSegments = Math.floor( radialSegments ) || 8;
	this.heightSegments = Math.floor( heightSegments ) || 1;
    revFn = revFn || function(x) { return 1; };

	openEnded = true;
	thetaStart = 0.0;
	thetaLength = Math.PI * 2;

	// buffers
	var indices = [];
	var vertices = [];
	var normals = [];
	var uvs = [];

	// helper variables
	var index = 0;
	var indexArray = [];
	var halfHeight = height / 2;
	var groupStart = 0;

	// generate geometry
	this.generateTorso = function(rfn) {
		var x, y;
		var normal = new THREE.Vector3();
		var vertex = new THREE.Vector3();

		this.groupCount = 0;

		// this will be used to calculate the normal
		var slope = 1;

		// generate vertices, normals and uvs

		for ( y = 0; y <= this.heightSegments; y ++ ) {

			var indexRow = [];

			var v = y / this.heightSegments;

			// calculate the radius of the current row

			var coords = rfn(v);
            var xx = coords[0], radius = coords[1];

			for ( x = 0; x <= this.radialSegments; x ++ ) {

				var u = x / this.radialSegments;

				var theta = u * thetaLength + thetaStart;

				var sinTheta = Math.sin( theta );
				var cosTheta = Math.cos( theta );

				// vertex

				vertex.x = radius * sinTheta;
				vertex.y = xx;
				vertex.z = radius * cosTheta;
				vertices.push( vertex.x, vertex.y, vertex.z );

				// normal

				normal.set( sinTheta, slope, cosTheta ).normalize();
				normals.push( normal.x, normal.y, normal.z );

				// uv

				uvs.push( u, 1 - v );

				// save index of vertex in respective row

				indexRow.push( index ++ );

			}

			// now save vertices of the row in our index array

			indexArray.push( indexRow );

		}

		// generate indices

		for ( x = 0; x < this.radialSegments; x ++ ) {

			for ( y = 0; y < this.heightSegments; y ++ ) {

				// we use the index array to access the correct indices

				var a = indexArray[ y ][ x ];
				var b = indexArray[ y + 1 ][ x ];
				var c = indexArray[ y + 1 ][ x + 1 ];
				var d = indexArray[ y ][ x + 1 ];

				// faces

				indices.push( a, b, d );
				indices.push( b, c, d );

				// update group counter

				this.groupCount += 6;

			}

		}

		// add a group to the geometry. this will ensure multi material support

		this.addGroup( this.groupStart, this.groupCount, 0 );

		// calculate new start value for groups

		this.groupStart += this.groupCount;

        this.setIndex( indices );
        this.addAttribute( 'position', new THREE.Float32BufferAttribute( vertices, 3 ) );
        this.addAttribute( 'normal', new THREE.Float32BufferAttribute( normals, 3 ) );
        this.addAttribute( 'uv', new THREE.Float32BufferAttribute( uvs, 2 ) );
	};
	this.generateTorso(revFn);
}

RevolutionBufferGeometry.prototype = Object.create( THREE.BufferGeometry.prototype );
RevolutionBufferGeometry.prototype.constructor = RevolutionBufferGeometry;
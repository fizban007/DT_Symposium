precision mediump float;

uniform sampler2D cubeTex;
uniform sampler2D transferTex;
uniform vec3 localCamPos;
uniform float steps;
uniform float alphaCorrection;
uniform float res;
uniform float row;

varying vec3 worldSpaceCoords;

// The maximum distance through our rendering volume is sqrt(3).
// The maximum number of steps we take to travel a distance of 1 is 512.
// ceil( sqrt(3) * 512 ) = 887
// This prevents the back of the image from getting cut off when steps=512 &
// viewing diagonally.
const int MAX_STEPS = 887;

// Acts like a texture3D using Z slices and trilinear filtering.
vec4 sampleAs3DTexture(vec3 texCoord) {
  vec4 colorSlice1;
  vec4 colorSlice2;
  vec2 texCoordSlice1;
  vec2 texCoordSlice2;

  // The z coordinate determines which Z slice we have to look for.
  // Z slice number goes from 0 to 255.
  float zSliceNumber1 = floor(texCoord.z * (res - 1.0));

  // As we use trilinear we go the next Z slice.
  float zSliceNumber2 = min(zSliceNumber1 + 1.0, (res - 1.0)); // Clamp to 255

  float col = res/row;
  // The Z slices are stored in a matrix of 16x16 of Z slices.
  // The original UV coordinates have to be rescaled by the tile numbers in each
  // row and column.
  texCoord.x /= row;
  texCoord.y /= col;

  texCoordSlice1 = texCoordSlice2 = texCoord.xy;

  // Add an offset to the original UV coordinates depending on the row and
  // column number.
  texCoordSlice1.x += (mod(zSliceNumber1, row) / row);
  texCoordSlice1.y += floor(zSliceNumber1 / row) / col;

  texCoordSlice2.x += (mod(zSliceNumber2, row) / row);
  texCoordSlice2.y += floor(zSliceNumber2 / row) / col;

  // Get the opacity value from the 2D texture.
  // Bilinear filtering is done at each texture2D by default.
  colorSlice1 = texture2D(cubeTex, texCoordSlice1);
  colorSlice2 = texture2D(cubeTex, texCoordSlice2);

  // TODO: only interpolate the sampled alpha values, then do the transfer lookup

  // Based on the opacity obtained earlier, get the RGB color in the transfer
  // function texture.
  // colorSlice1.rgb = texture2D(transferTex, vec2(colorSlice1.a, 1.0)).rgb;
  // colorSlice2.rgb = texture2D(transferTex, vec2(colorSlice2.a, 1.0)).rgb;
  colorSlice1.rgb = texture2D(transferTex, vec2(colorSlice1.a, 2.0)).rgb;
  colorSlice2.rgb = texture2D(transferTex, vec2(colorSlice2.a, 2.0)).rgb;
  // float c1 = (colorSlice1.r - colorSlice1.b) / 0.5 + 0.5;
  // float c2 = (colorSlice2.r - colorSlice2.b) / 0.5 + 0.5;
  // colorSlice1.rgb = texture2D(transferTex, vec2(c1, 1.0)).rgb;
  // colorSlice2.rgb = texture2D(transferTex, vec2(c2, 1.0)).rgb;
  // colorSlice1.rgb = vec4(colorSlice1.r, 0.01, colorSlice1.b, 1.0).rgb;
  // colorSlice2.rgb = vec4(colorSlice2.r, 0.01, colorSlice2.b, 1.0).rgb;

  // How distant is zSlice1 to ZSlice2. Used to interpolate between one Z slice
  // and the other.
  float zDifference = mod(texCoord.z * (res - 1.0), 1.0);

  // Finally interpolate between the two intermediate colors of each Z slice.
  return mix(colorSlice1, colorSlice2, zDifference);
}

struct AABB {
    vec3 min; // (0,0,0) 
    vec3 max; // (1,1,1)
};

struct Ray {
    vec3 origin; vec3 dir;
};

float getUnitAABBEntry( in Ray r ) {
   AABB b;
   b.min = vec3( 0 ); 
   b.max = vec3( 1 );

   // compute clipping for box.min and box.max corner
   vec3 rInvDir = vec3( 1.0 ) / r.dir;
   vec3 tMinima = ( b.min - r.origin ) * rInvDir; 
   vec3 tMaxima = ( b.max - r.origin ) * rInvDir;

   // sort for nearest corner
   vec3 tEntries = min( tMinima, tMaxima );

   // find first real entry value of 3 t-distance values in vec3 container
   vec2 tMaxEntryCandidates = max( vec2( tEntries.st ), vec2( tEntries.pp ) ); 
   // float tMaxEntry = max( tMaxEntryCandidates.s, tMaxEntryCandidates.t );
   return max( tMaxEntryCandidates.s, tMaxEntryCandidates.t );
}

vec3 getCloserPos( in vec3 camera, in vec3 frontFaceIntersection, in float t ) {
    float useFrontCoord = 0.5 + 0.5 * sign( t );
    vec3 startPos = mix( camera, frontFaceIntersection, useFrontCoord );   
    return startPos;
}

void main(void)
{
    Ray r;
    r.origin = localCamPos;
    r.dir = normalize( worldSpaceCoords - localCamPos );

    float t = getUnitAABBEntry( r );
    vec3 frontFaceLocalUnitTexCoord = r.origin + r.dir * t;
    vec3 startPos = getCloserPos( localCamPos, frontFaceLocalUnitTexCoord, t );

    // loop for integration follows here
    vec3 currentPosition = startPos;
    vec3 end = worldSpaceCoords;
    float rayLength = length(end - currentPosition);
    float delta = 1.0 / steps;

    // The color accumulator.
    vec4 accumulatedColor = vec4(0.0);

    // The alpha value accumulated so far.
    float accumulatedAlpha = 0.0;

    // How long has the ray travelled so far.
    float accumulatedLength = 0.0;

    // If we have twice as many samples, we only need ~1/2 the alpha per sample.
    // Scaling by 256/10 just happens to give a good value for the alphaCorrection
    // slider.
    float alphaScaleFactor = 25.6 * delta;

    vec4 colorSample;
    float alphaSample;
    
    // Perform the ray marching iterations
    for (int i = 0; i < MAX_STEPS; i++) {
      // Get the voxel intensity value from the 3D texture.
      colorSample = sampleAs3DTexture(currentPosition);

      // Allow the alpha correction customization.
      alphaSample = colorSample.a * alphaCorrection;

      // Applying this effect to both the color and alpha accumulation results in
      // more realistic transparency.
      alphaSample *= (1.0 - accumulatedAlpha);

      // Scaling alpha by the number of steps makes the final color invariant to
      // the step size.
      alphaSample *= alphaScaleFactor;

      // Perform the composition.
      accumulatedColor += colorSample * alphaSample;

      // Store the alpha accumulated so far.
      accumulatedAlpha += alphaSample;

      // Advance the ray.
      currentPosition += r.dir * delta;
      accumulatedLength += delta;

      // If the length traversed is more than the ray length, or if the alpha
      // accumulated reaches 1.0 then exit.
      if (accumulatedLength >= rayLength || accumulatedAlpha >= 1.0)
        break;
    }

    gl_FragColor = accumulatedColor;
    // gl_FragColor = vec4(0.0, 0.2, 0.2, 1.0);

}

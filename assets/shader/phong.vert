#version 330 core
layout(location = 0) in vec3 Position_in;
layout(location = 1) in vec3 Normal_in;
layout(location = 2) in vec2 TextureCoordinate_in;

out VS_OUT {
  vec3 Position;
  vec3 Normal;
  vec2 TextureCoordinate;
  vec4 LightSpacePosition;
  flat vec4 lightVector;
  flat vec4 viewPosition;
  flat vec4 lightCoefficients;
  vec3 rawPosition;
} vs_out;

layout (std140) uniform model {
  // Model matrix
  mat4 modelMatrix;
  // inverse(transpose(model)), precalculate using CPU for efficiency
  mat4 normalMatrix;
};

layout (std140) uniform camera {
  // Projection * View matrix
  mat4 viewProjectionMatrix;
  // Position of the camera
  vec4 viewPosition;
};

layout (std140) uniform light {
  // Projection * View matrix
  mat4 lightSpaceMatrix;
  // Position or direction of the light
  vec4 lightVector;
  // inner cutoff, outer cutoff, isSpotlight, isDirectionalLight
  vec4 coefficients;
};

void main() {
  // TODO: vertex shader / fragment shader
  // Hint:
  //       1. how to write a vertex shader:
  //          a. The output is gl_Position and anything you want to pass to the fragment shader. (Apply matrix multiplication yourself)
  //       2. how to write a fragment shader:
  //          a. The output is FragColor (any var is OK)
  //       3. colors
  //          a. For point light & directional light, lighting = ambient + attenuation * shadow * (diffuse + specular)
  //          b. If you want to implement multiple light sources, you may want to use lighting = shadow * attenuation * (ambient + (diffuse + specular))
  //       4. attenuation
  //          a. spotlight & pointlight: see spec
  //          b. directional light = no
  //          c. Use formula from slides 'shading.ppt' page 20
  //       5. spotlight cutoff: inner and outer from coefficients.x and coefficients.y
  //       6. diffuse = kd * max(normal vector dot light direction, 0.0)
  //       7. specular = ks * pow(max(normal vector dot halfway direction), 0.0), 8.0);
  //       8. notice the difference of light direction & distance between directional light & point light
  //       9. we've set ambient & color for you
  vs_out.Position = vec3(modelMatrix * vec4(Position_in, 1.0));
  vs_out.Normal = mat3(normalMatrix) * Normal_in;
  vs_out.TextureCoordinate = TextureCoordinate_in;
  vs_out.LightSpacePosition = lightSpaceMatrix * vec4(vs_out.Position, 1.0);
  vs_out.lightVector = lightVector;
  vs_out.viewPosition = viewPosition;
  vs_out.lightCoefficients = coefficients;
  vs_out.rawPosition = mat3(modelMatrix) * Position_in;
  gl_Position = viewProjectionMatrix * modelMatrix * vec4(Position_in, 1.0);
}

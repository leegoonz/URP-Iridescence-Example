#ifndef UNIVERSAL_BRDF_INCLUDED
#define UNIVERSAL_BRDF_INCLUDED

#include "Packages/com.unity.render-pipelines.core/ShaderLibrary/BSDF.hlsl"
#include "Packages/com.unity.render-pipelines.core/ShaderLibrary/CommonMaterial.hlsl"
#include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/Deprecated.hlsl"
#include "Packages/com.unity.render-pipelines.universal/ShaderLibrary/SurfaceData.hlsl"

#define kDielectricSpec half4(0.04, 0.04, 0.04, 1.0 - 0.04) // standard dielectric reflectivity coef at incident angle (= 4%)

///////////////////////////////////////////////////////////////////////////////
//                         Helper Functions                                  //
///////////////////////////////////////////////////////////////////////////////

// Common constants
#define PI 3.14159265358979323846

// XYZ to CIE 1931 RGB color space (using neutral E illuminant)
static const half3x3 XYZ_TO_RGB = half3x3(2.3706743, -0.5138850, 0.0052982, -0.9000405, 1.4253036, -0.0146949, -0.4706338, 0.0885814, 1.0093968);

// Square functions for cleaner code
inline float sqr(float x) { return x * x; }
inline float2 sqr(float2 x) { return x * x; }

// Depolarization functions for natural light
inline float depol(float2 polV) { return 0.5 * (polV.x + polV.y); }
inline float3 depolColor(float3 colS, float3 colP) { return 0.5 * (colS + colP); }

// Evaluation XYZ sensitivity curves in Fourier space
float3 evalSensitivity(float opd, float shift) {

    // Use Gaussian fits, given by 3 parameters: val, pos and var
    float phase = 2 * PI * opd * 1.0e-6;
    float3 val = float3(5.4856e-13, 4.4201e-13, 5.2481e-13);
    float3 pos = float3(1.6810e+06, 1.7953e+06, 2.2084e+06);
    float3 var = float3(4.3278e+09, 9.3046e+09, 6.6121e+09);
    float3 xyz = val * sqrt(2.0 * PI * var) * cos(pos * phase + shift) * exp(-var * phase * phase);
    xyz.x += 9.7470e-14 * sqrt(2.0 * PI * 4.5282e+09) * cos(2.2399e+06 * phase + shift) * exp(-4.5282e+09 * phase * phase);
    return xyz / 1.0685e-7;
}

///////////////////////////////////////////////////////////////////////////////
//                         BRDF Functions                                    //
///////////////////////////////////////////////////////////////////////////////

struct BRDFData
{
    half3 albedo;
    half3 diffuse;
    half3 specular;
    half reflectivity;
    half perceptualRoughness;
    half roughness;
    half roughness2;
    half grazingTerm;

    // We save some light invariant BRDF terms so we don't have to recompute
    // them in the light loop. Take a look at DirectBRDF function for detailed explaination.
    half normalizationTerm;     // roughness * 4.0 + 2.0
    half roughness2MinusOne;    // roughness^2 - 1.0
#ifdef _IRIDESCENCE
    half iridescenceThickness;
    half iridescenceEta2;
    half iridescenceEta3;
    half iridescenceKappa3;
#endif
};

half ReflectivitySpecular(half3 specular)
{
#if defined(SHADER_API_GLES)
    return specular.r; // Red channel - because most metals are either monocrhome or with redish/yellowish tint
#else
    return Max3(specular.r, specular.g, specular.b);
#endif
}

half OneMinusReflectivityMetallic(half metallic)
{
    // We'll need oneMinusReflectivity, so
    //   1-reflectivity = 1-lerp(dielectricSpec, 1, metallic) = lerp(1-dielectricSpec, 0, metallic)
    // store (1-dielectricSpec) in kDielectricSpec.a, then
    //   1-reflectivity = lerp(alpha, 0, metallic) = alpha + metallic*(0 - alpha) =
    //                  = alpha - metallic * alpha
    half oneMinusDielectricSpec = kDielectricSpec.a;
    return oneMinusDielectricSpec - metallic * oneMinusDielectricSpec;
}

half MetallicFromReflectivity(half reflectivity)
{
    half oneMinusDielectricSpec = kDielectricSpec.a;
    return (reflectivity - kDielectricSpec.r) / oneMinusDielectricSpec;
}

inline void InitializeBRDFDataDirect(half3 albedo, half3 diffuse, half3 specular, half reflectivity, half oneMinusReflectivity, half smoothness, inout half alpha,
#ifdef _IRIDESCENCE
    half iridescenceThickness,
    half iridescenceEta2,
    half iridescenceEta3,
    half iridescenceKappa3,
#endif
out BRDFData outBRDFData)
{
    outBRDFData = (BRDFData)0;
    outBRDFData.albedo = albedo;
    outBRDFData.diffuse = diffuse;
    outBRDFData.specular = specular;
    outBRDFData.reflectivity = reflectivity;

    outBRDFData.perceptualRoughness = PerceptualSmoothnessToPerceptualRoughness(smoothness);
    outBRDFData.roughness           = max(PerceptualRoughnessToRoughness(outBRDFData.perceptualRoughness), HALF_MIN_SQRT);
    outBRDFData.roughness2          = max(outBRDFData.roughness * outBRDFData.roughness, HALF_MIN);
    outBRDFData.grazingTerm         = saturate(smoothness + reflectivity);
    outBRDFData.normalizationTerm   = outBRDFData.roughness * half(4.0) + half(2.0);
    outBRDFData.roughness2MinusOne  = outBRDFData.roughness2 - half(1.0);

#ifdef _IRIDESCENCE
    outBRDFData.iridescenceThickness = iridescenceThickness;
    outBRDFData.iridescenceEta2 = iridescenceEta2;
    outBRDFData.iridescenceEta3 = iridescenceEta3;
    outBRDFData.iridescenceKappa3 = iridescenceKappa3;
#endif
    

#ifdef _ALPHAPREMULTIPLY_ON
    outBRDFData.diffuse *= alpha;
    alpha = alpha * oneMinusReflectivity + reflectivity; // NOTE: alpha modified and propagated up.
#endif
}

// Legacy: do not call, will not correctly initialize albedo property.
inline void InitializeBRDFDataDirect(half3 diffuse, half3 specular, half reflectivity, half oneMinusReflectivity, half smoothness, inout half alpha,
#ifdef _IRIDESCENCE
    half iridescenceThickness,
    half iridescenceEta2,
    half iridescenceEta3,
    half iridescenceKappa3,
#endif
out BRDFData outBRDFData)
{
    InitializeBRDFDataDirect(half3(0.0, 0.0, 0.0), diffuse, specular, reflectivity, oneMinusReflectivity, smoothness, alpha,
    #ifdef _IRIDESCENCE
        iridescenceThickness,
        iridescenceEta2,
        iridescenceEta3,
        iridescenceKappa3,
    #endif
    outBRDFData);
}

// Initialize BRDFData for material, managing both specular and metallic setup using shader keyword _SPECULAR_SETUP.
inline void InitializeBRDFData(half3 albedo, half metallic, half3 specular, half smoothness, inout half alpha,
#ifdef _IRIDESCENCE
    half iridescenceThickness,
    half iridescenceEta2,
    half iridescenceEta3,
    half iridescenceKappa3,
#endif
out BRDFData outBRDFData)
{
#ifdef _SPECULAR_SETUP
    half reflectivity = ReflectivitySpecular(specular);
    half oneMinusReflectivity = half(1.0) - reflectivity;
    half3 brdfDiffuse = albedo * (half3(1.0, 1.0, 1.0) - specular);
    half3 brdfSpecular = specular;
#else
    half oneMinusReflectivity = OneMinusReflectivityMetallic(metallic);
    half reflectivity = half(1.0) - oneMinusReflectivity;
    half3 brdfDiffuse = albedo * oneMinusReflectivity;
    half3 brdfSpecular = lerp(kDieletricSpec.rgb, albedo, metallic);
#endif

    InitializeBRDFDataDirect(albedo, brdfDiffuse, brdfSpecular, reflectivity, oneMinusReflectivity, smoothness, alpha,
    #ifdef _IRIDESCENCE
        iridescenceThickness,
        iridescenceEta2,
        iridescenceEta3,
        iridescenceKappa3,
    #endif
    outBRDFData);
}

inline void InitializeBRDFData(inout SurfaceData surfaceData, out BRDFData brdfData)
{
    InitializeBRDFData(surfaceData.albedo, surfaceData.metallic, surfaceData.specular, surfaceData.smoothness, surfaceData.alpha,
#ifdef _IRIDESCENCE
surfaceData.iridescenceThickness,
surfaceData.iridescenceEta2,
surfaceData.iridescenceEta3,
surfaceData.iridescenceKappa3,
#endif
    brdfData);
}

half3 ConvertF0ForClearCoat15(half3 f0)
{
#if defined(SHADER_API_MOBILE)
    return ConvertF0ForAirInterfaceToF0ForClearCoat15Fast(f0);
#else
    return ConvertF0ForAirInterfaceToF0ForClearCoat15(f0);
#endif
}


inline void InitializeBRDFDataClearCoat(half clearCoatMask, half clearCoatSmoothness, inout BRDFData baseBRDFData, out BRDFData outBRDFData)
{
    outBRDFData = (BRDFData)0;
    outBRDFData.albedo = half(1.0);

    // Calculate Roughness of Clear Coat layer
    outBRDFData.diffuse             = kDielectricSpec.aaa; // 1 - kDielectricSpec
    outBRDFData.specular            = kDielectricSpec.rgb;
    outBRDFData.reflectivity        = kDielectricSpec.r;

    outBRDFData.perceptualRoughness = PerceptualSmoothnessToPerceptualRoughness(clearCoatSmoothness);
    outBRDFData.roughness           = max(PerceptualRoughnessToRoughness(outBRDFData.perceptualRoughness), HALF_MIN_SQRT);
    outBRDFData.roughness2          = max(outBRDFData.roughness * outBRDFData.roughness, HALF_MIN);
    outBRDFData.normalizationTerm   = outBRDFData.roughness * half(4.0) + half(2.0);
    outBRDFData.roughness2MinusOne  = outBRDFData.roughness2 - half(1.0);
    outBRDFData.grazingTerm         = saturate(clearCoatSmoothness + kDielectricSpec.x);

// Relatively small effect, cut it for lower quality
#if !defined(SHADER_API_MOBILE)
    // Modify Roughness of base layer using coat IOR
    half ieta                        = lerp(1.0h, CLEAR_COAT_IETA, clearCoatMask);
    half coatRoughnessScale          = Sq(ieta);
    half sigma                       = RoughnessToVariance(PerceptualRoughnessToRoughness(baseBRDFData.perceptualRoughness));

    baseBRDFData.perceptualRoughness = RoughnessToPerceptualRoughness(VarianceToRoughness(sigma * coatRoughnessScale));

    // Recompute base material for new roughness, previous computation should be eliminated by the compiler (as it's unused)
    baseBRDFData.roughness          = max(PerceptualRoughnessToRoughness(baseBRDFData.perceptualRoughness), HALF_MIN_SQRT);
    baseBRDFData.roughness2         = max(baseBRDFData.roughness * baseBRDFData.roughness, HALF_MIN);
    baseBRDFData.normalizationTerm  = baseBRDFData.roughness * 4.0h + 2.0h;
    baseBRDFData.roughness2MinusOne = baseBRDFData.roughness2 - 1.0h;
#endif

    // Darken/saturate base layer using coat to surface reflectance (vs. air to surface)
    baseBRDFData.specular = lerp(baseBRDFData.specular, ConvertF0ForClearCoat15(baseBRDFData.specular), clearCoatMask);
    // TODO: what about diffuse? at least in specular workflow diffuse should be recalculated as it directly depends on it.
}

BRDFData CreateClearCoatBRDFData(SurfaceData surfaceData, inout BRDFData brdfData)
{
    BRDFData brdfDataClearCoat = (BRDFData)0;

    #if defined(_CLEARCOAT) || defined(_CLEARCOATMAP)
    // base brdfData is modified here, rely on the compiler to eliminate dead computation by InitializeBRDFData()
    InitializeBRDFDataClearCoat(surfaceData.clearCoatMask, surfaceData.clearCoatSmoothness, brdfData, brdfDataClearCoat);
    #endif
    return brdfDataClearCoat;
}

// Computes the specular term for EnvironmentBRDF
half3 EnvironmentBRDFSpecular(BRDFData brdfData, half fresnelTerm)
{
    float surfaceReduction = 1.0 / (brdfData.roughness2 + 1.0);
    return half3(surfaceReduction * lerp(brdfData.specular, brdfData.grazingTerm, fresnelTerm));
}

half3 EnvironmentBRDF(BRDFData brdfData, half3 indirectDiffuse, half3 indirectSpecular, half fresnelTerm)
{
    half3 c = indirectDiffuse * brdfData.diffuse;
    c += indirectSpecular * EnvironmentBRDFSpecular(brdfData, fresnelTerm);
    return c;
}

half3 EnvironmentBRDFIridescence(BRDFData brdfData, half3 indirectDiffuse, half3 indirectSpecular, half3 fresnelIridescent)
{
    half3 c = indirectDiffuse * brdfData.diffuse;
    float surfaceReduction = 1.0 / (brdfData.roughness2 + 1.0);
    c += surfaceReduction * indirectSpecular * lerp(brdfData.specular * fresnelIridescent, brdfData.grazingTerm, fresnelIridescent);
    return c;
}

// Environment BRDF without diffuse for clear coat
half3 EnvironmentBRDFClearCoat(BRDFData brdfData, half clearCoatMask, half3 indirectSpecular, half fresnelTerm)
{
    float surfaceReduction = 1.0 / (brdfData.roughness2 + 1.0);
    return indirectSpecular * EnvironmentBRDFSpecular(brdfData, fresnelTerm) * clearCoatMask;
}


#ifdef _IRIDESCENCE
// GGX distribution function
float GGX(float NdotH, float a) {
    float a2 = sqr(a);
    return a2 / (PI * sqr(sqr(NdotH) * (a2 - 1) + 1));
}

// Smith GGX geometric functions
float smithG1_GGX(float NdotV, float a) {
    float a2 = sqr(a);
    return 2 / (1 + sqrt(1 + a2 * (1 - sqr(NdotV)) / sqr(NdotV)));
}

float smithG_GGX(float NdotL, float NdotV, float a) {
    return smithG1_GGX(NdotL, a) * smithG1_GGX(NdotV, a);
}

// Fresnel equations for dielectric/dielectric interfaces.
void fresnelDielectric(in float ct1, in float n1, in float n2,
    out float2 R, out float2 phi) {

    float st1 = (1 - ct1 * ct1); // Sinus theta1 'squared'
    float nr = n1 / n2;

    if (sqr(nr) * st1 > 1) { // Total reflection

        float2 R = float2(1, 1);
        phi = 2.0 * atan(float2(-sqr(nr) * sqrt(st1 - 1.0 / sqr(nr)) / ct1,
            -sqrt(st1 - 1.0 / sqr(nr)) / ct1));
    }
    else {   // Transmission & Reflection

        float ct2 = sqrt(1 - sqr(nr) * st1);
        float2 r = float2((n2 * ct1 - n1 * ct2) / (n2 * ct1 + n1 * ct2),
            (n1 * ct1 - n2 * ct2) / (n1 * ct1 + n2 * ct2));
        phi.x = (r.x < 0.0) ? PI : 0.0;
        phi.y = (r.y < 0.0) ? PI : 0.0;
        R = sqr(r);
    }
}

// Fresnel equations for dielectric/conductor interfaces.
void fresnelConductor(in float ct1, in float n1, in float n2, in float k,
    out float2 R, out float2 phi) {

    if (k == 0) { // use dielectric formula to avoid numerical issues
        fresnelDielectric(ct1, n1, n2, R, phi);
        return;
    }

    float A = sqr(n2) * (1 - sqr(k)) - sqr(n1) * (1 - sqr(ct1));
    float B = sqrt(sqr(A) + sqr(2 * sqr(n2) * k));
    float U = sqrt((A + B) / 2.0);
    float V = sqrt((B - A) / 2.0);

    R.y = (sqr(n1 * ct1 - U) + sqr(V)) / (sqr(n1 * ct1 + U) + sqr(V));
    phi.y = atan2(2 * n1 * V * ct1, sqr(U) + sqr(V) - sqr(n1 * ct1)) + PI;

    R.x = (sqr(sqr(n2) * (1 - sqr(k)) * ct1 - n1 * U) + sqr(2 * sqr(n2) * k * ct1 - n1 * V))
        / (sqr(sqr(n2) * (1 - sqr(k)) * ct1 + n1 * U) + sqr(2 * sqr(n2) * k * ct1 + n1 * V));
    phi.x = atan2(2 * n1 * sqr(n2) * ct1 * (2 * k * U - (1 - sqr(k)) * V), sqr(sqr(n2) * (1 + sqr(k)) * ct1) - sqr(n1) * (sqr(U) + sqr(V)));
}


// Evaluate the reflectance for a thin-film layer on top of a dielectric medum
// Based on the paper [LAURENT 2017] A Practical Extension to Microfacet Theory for the Modeling of Varying Iridescence
half3 ThinFilmIridescence(BRDFData brdfData, float cosTheta1)
{
    float eta_1 = 1.0; // Air on top, no coat.
    float eta_2 = brdfData.iridescenceEta2;
    float eta_3 = brdfData.iridescenceEta3;
    float kappa_3 = brdfData.iridescenceKappa3;

    // iridescenceThickness unit is micrometer for this equation here. Mean 0.5 is 500nm.
    float Dinc = 2 * eta_2 * brdfData.iridescenceThickness;

    // Force eta_2 -> eta_1 when Dinc -> 0.0
    eta_2 = lerp(eta_1, eta_2, smoothstep(0.0, 0.03, Dinc));

    float cosTheta2 = sqrt(1.0 - sqr(eta_1 / eta_2) * (1 - sqr(cosTheta1)));

    // First interface
    float2 R12, phi12;
    fresnelDielectric(cosTheta1, eta_1, eta_2, R12, phi12);
    float2 R21 = R12;
    float2 T121 = float2(1.0, 1.0) - R12;
    float2 phi21 = float2(PI, PI) - phi12;

    // Second interface
    float2 R23, phi23;
    fresnelConductor(cosTheta2, eta_2, eta_3, kappa_3, R23, phi23);

    // Phase shift
    float OPD = Dinc * cosTheta2;
    float2 phi2 = phi21 + phi23;

    // Compound terms
    float3 I = float3(0, 0, 0);
    float2 R123 = clamp(R12 * R23, 1e-5, 0.9999);
    float2 r123 = sqrt(R123);
    float2 Rs = sqr(T121) * R23 / (float2(1.0, 1.0) - R123);

    // Reflectance term for m=0 (DC term amplitude)
    float2 C0 = R12 + Rs;
    float3 S0 = evalSensitivity(0.0, 0.0);
    I += depol(C0) * S0;

    // Reflectance term for m>0 (pairs of diracs)
    float2 Cm = Rs - T121;

    [unroll(3)]
    for (int m = 1; m <= 3; ++m)
    {
        Cm *= r123;
        float3 SmS = 2.0 * evalSensitivity(m * OPD, m * phi2.x);
        float3 SmP = 2.0 * evalSensitivity(m * OPD, m * phi2.y);
        I += depolColor(Cm.x * SmS, Cm.y * SmP);
    }

    // Convert back to RGB reflectance
    I = max(mul(I, XYZ_TO_RGB), float3(0.0, 0.0, 0.0));

    return I;
}

half3 DirectBDRFIridescence(BRDFData brdfData, half3 normalWS, half3 lightDirectionWS, half3 viewDirectionWS)
{
    // Compute dot products
    float NdotL = dot(normalWS, lightDirectionWS);
    float NdotV = dot(normalWS, viewDirectionWS);

    float3 halfDir = SafeNormalize(float3(lightDirectionWS) + float3(viewDirectionWS));
    float NdotH = dot(normalWS, halfDir);
    float cosTheta1 = dot(halfDir, float3(lightDirectionWS));

    half3 I = ThinFilmIridescence(brdfData, cosTheta1);
    // Microfacet BRDF formula
    float D = GGX(NdotH, brdfData.perceptualRoughness);
    float G = smithG_GGX(NdotL, NdotV, brdfData.perceptualRoughness);
    half3 specularTerm = D * G * I / (4 * NdotL * NdotV);
    return specularTerm;
}
#endif

// Computes the scalar specular term for Minimalist CookTorrance BRDF
// NOTE: needs to be multiplied with reflectance f0, i.e. specular color to complete
half DirectBRDFSpecular(BRDFData brdfData, half3 normalWS, half3 lightDirectionWS, half3 viewDirectionWS)
{
    float3 lightDirectionWSFloat3 = float3(lightDirectionWS);
    float3 halfDir = SafeNormalize(lightDirectionWSFloat3 + float3(viewDirectionWS));

    float NoH = saturate(dot(float3(normalWS), halfDir));
    half LoH = half(saturate(dot(lightDirectionWSFloat3, halfDir)));

    // GGX Distribution multiplied by combined approximation of Visibility and Fresnel
    // BRDFspec = (D * V * F) / 4.0
    // D = roughness^2 / ( NoH^2 * (roughness^2 - 1) + 1 )^2
    // V * F = 1.0 / ( LoH^2 * (roughness + 0.5) )
    // See "Optimizing PBR for Mobile" from Siggraph 2015 moving mobile graphics course
    // https://community.arm.com/events/1155

    // Final BRDFspec = roughness^2 / ( NoH^2 * (roughness^2 - 1) + 1 )^2 * (LoH^2 * (roughness + 0.5) * 4.0)
    // We further optimize a few light invariant terms
    // brdfData.normalizationTerm = (roughness + 0.5) * 4.0 rewritten as roughness * 4.0 + 2.0 to a fit a MAD.
    float d = NoH * NoH * brdfData.roughness2MinusOne + 1.00001f;
    half d2 = half(d * d);

    half LoH2 = LoH * LoH;
    half specularTerm = brdfData.roughness2 / (d2 * max(half(0.1), LoH2) * brdfData.normalizationTerm);

    // On platforms where half actually means something, the denominator has a risk of overflow
    // clamp below was added specifically to "fix" that, but dx compiler (we convert bytecode to metal/gles)
    // sees that specularTerm have only non-negative terms, so it skips max(0,..) in clamp (leaving only min(100,...))
#if defined (SHADER_API_MOBILE) || defined (SHADER_API_SWITCH)
    specularTerm = specularTerm - HALF_MIN;
    specularTerm = clamp(specularTerm, 0.0, 100.0); // Prevent FP16 overflow on mobiles
#endif

return specularTerm;
}

// Based on Minimalist CookTorrance BRDF
// Implementation is slightly different from original derivation: http://www.thetenthplanet.de/archives/255
//
// * NDF [Modified] GGX
// * Modified Kelemen and Szirmay-Kalos for Visibility term
// * Fresnel approximated with 1/LdotH
half3 DirectBDRF(BRDFData brdfData, half3 normalWS, half3 lightDirectionWS, half3 viewDirectionWS, bool specularHighlightsOff)
{
    // Can still do compile-time optimisation.
    // If no compile-time optimized, extra overhead if branch taken is around +2.5% on some untethered platforms, -10% if not taken.
    [branch] if (!specularHighlightsOff)
    {
        half specularTerm = 0;
#ifdef _IRIDESCENCE
        specularTerm = DirectBDRFIridescence(brdfData, normalWS, lightDirectionWS, viewDirectionWS);
#else
        specularTerm = DirectBRDFSpecular(brdfData, normalWS, lightDirectionWS, viewDirectionWS);
#endif
        half3 color = brdfData.diffuse + specularTerm * brdfData.specular;
        return color;
    }
    else
        return brdfData.diffuse;
}

// Based on Minimalist CookTorrance BRDF
// Implementation is slightly different from original derivation: http://www.thetenthplanet.de/archives/255
//
// * NDF [Modified] GGX
// * Modified Kelemen and Szirmay-Kalos for Visibility term
// * Fresnel approximated with 1/LdotH
half3 DirectBRDF(BRDFData brdfData, half3 normalWS, half3 lightDirectionWS, half3 viewDirectionWS)
{
#ifndef _SPECULARHIGHLIGHTS_OFF
    return brdfData.diffuse + DirectBRDFSpecular(brdfData, normalWS, lightDirectionWS, viewDirectionWS) * brdfData.specular;
#else
    return brdfData.diffuse;
#endif
}

#endif

using System;

namespace UnityEngine.Rendering.Universal
{
    public enum TonemappingMode
    {
        None,
        Neutral, // Neutral tonemapper
        ACES,    // ACES Filmic reference tonemapper (custom approximation)
        XD_ACES,//XRP_FEATURE ACES UE4 Log Adding here
    }

    [Serializable, VolumeComponentMenu("Post-processing/Tonemapping")]
    public sealed class Tonemapping : VolumeComponent, IPostProcessComponent
    {
        [Tooltip("Select a tonemapping algorithm to use for the color grading process.")]
        public TonemappingModeParameter mode = new TonemappingModeParameter(TonemappingMode.None);

        //XRP_FEATURE
        //create floatParameter
        [Tooltip("Film_Slope")]
        public ClampedFloatParameter slope = new ClampedFloatParameter(0.88f,0f,1f);
        [Tooltip("Film_Toe")]
        public ClampedFloatParameter toe = new ClampedFloatParameter(0.55f, 0f, 1f);
        [Tooltip("Film_Shoulder")]
        public ClampedFloatParameter shoulder = new ClampedFloatParameter(0.26f, 0f, 1f);
        [Tooltip("Film_BlackClip")]
        public ClampedFloatParameter blackClip = new ClampedFloatParameter(0.0f, 0f, 1f);
        [Tooltip("Film_WhiteClip")]
        public ClampedFloatParameter whiteClip = new ClampedFloatParameter(0.04f, 0f, 1f);
        //XRP_FEATURE END

        public bool IsActive() => mode.value != TonemappingMode.None;

        public bool IsTileCompatible() => true;
    }

    [Serializable]
    public sealed class TonemappingModeParameter : VolumeParameter<TonemappingMode> { public TonemappingModeParameter(TonemappingMode value, bool overrideState = false) : base(value, overrideState) { } }
}

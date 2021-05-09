using UnityEngine.Rendering.Universal;
using UnityEngine;
using UnityEngine.Rendering;
namespace UnityEditor.Rendering.Universal
{
    [VolumeComponentEditor(typeof(Tonemapping))]
    sealed class TonemappingEditor : VolumeComponentEditor
    {
        SerializedDataParameter m_Mode;

        //XRP_FEATURE
        SerializedDataParameter m_Slope;
        SerializedDataParameter m_Toe;
        SerializedDataParameter m_Shoulder;
        SerializedDataParameter m_BlackClip;
        SerializedDataParameter m_WhiteClip;
        Tonemapping component;
        //XRP_FEATURE End



        public override void OnEnable()
        {
        //XRP_FEATURE
            component = (Tonemapping)target;
        //XRP_FEATURE End
        
            var o = new PropertyFetcher<Tonemapping>(serializedObject);
            m_Mode = Unpack(o.Find(x => x.mode));

        //XRP_FEATURE
            m_Slope      = Unpack(o.Find(x => x.slope));
            m_Toe        = Unpack(o.Find(x => x.toe));
            m_Shoulder   = Unpack(o.Find(x => x.shoulder));
            m_BlackClip  = Unpack(o.Find(x => x.blackClip));
            m_WhiteClip  = Unpack(o.Find(x => x.whiteClip));
        //XRP_FEATURE End


        }

        public override void OnInspectorGUI()
        {
            if (UniversalRenderPipeline.asset?.postProcessingFeatureSet == PostProcessingFeatureSet.PostProcessingV2)
            {
                EditorGUILayout.HelpBox(UniversalRenderPipelineAssetEditor.Styles.postProcessingGlobalWarning, MessageType.Warning);
                return;
            }
            PropertyField(m_Mode);

        //XRP_FEATURE
            if (component.mode.overrideState && component.mode == TonemappingMode.XD_ACES)
            {

                GUILayout.BeginHorizontal();
                PropertyField(m_Slope);
                if (GUILayout.Button("Reset"))
                {
                    component.slope.value = 0.88f;
                }
                GUILayout.EndHorizontal();

                GUILayout.BeginHorizontal();
                PropertyField(m_Toe);
                if (GUILayout.Button("Reset"))
                {
                    component.toe.value = 0.55f;
                }
                GUILayout.EndHorizontal();

                GUILayout.BeginHorizontal();
                PropertyField(m_Shoulder);
                if (GUILayout.Button("Reset"))
                {
                    component.shoulder.value = 0.26f;
                }
                GUILayout.EndHorizontal();

                GUILayout.BeginHorizontal();
                PropertyField(m_BlackClip);
                if (GUILayout.Button("Reset"))
                {
                    component.blackClip.value = 0.0f;
                }
                GUILayout.EndHorizontal();

                GUILayout.BeginHorizontal();
                PropertyField(m_WhiteClip);
                if (GUILayout.Button("Reset"))
                {
                    component.whiteClip.value = 0.04f;
                }
                GUILayout.EndHorizontal();
            }
        //XRP_FEATURE End


            // Display a warning if the user is trying to use a tonemap while rendering in LDR
            if (UniversalRenderPipeline.asset?.supportsHDR == false)
                EditorGUILayout.HelpBox("Tonemapping should only be used when working in HDR.", MessageType.Warning);
        }
    }
}

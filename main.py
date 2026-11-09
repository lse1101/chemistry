import streamlit as st
import requests
import pandas as pd
import plotly.express as px
from scipy.stats import pearsonr

# 페이지 기본 설정
st.set_page_config(
    page_title="탄소 사슬 길이에 따른 알코올 물성 변화 분석",
    page_icon="🧪",
    layout="wide"
)

st.title("🧪 탄소 사슬 길이에 따른 알코올의 분자 물성 변화 분석")
st.caption("PubChem PUG-REST API를 활용한 데이터 수집 및 시각화")
st.markdown("---")

# 분석 대상 데이터 정의 (탄소 수 기준 정렬)
ALCOHOLS = [
    {"name": "Methanol", "kr_name": "메탄올", "carbons": 1},
    {"name": "Ethanol", "kr_name": "에탄올", "carbons": 2},
    {"name": "1-Propanol", "kr_name": "1-프로판올", "carbons": 3},
    {"name": "1-Butanol", "kr_name": "1-부탄올", "carbons": 4},
    {"name": "1-Pentanol", "kr_name": "1-펜탄올", "carbons": 5},
]

# PubChem API 데이터 수집 함수 (캐싱 적용)
@st.cache_data(ttl=3600)
def fetch_pubchem_data(compound_name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/property/MolecularWeight,XLogP/JSON"
    try:
        response = requests.get(url, timeout=5)
        if response.status_code == 200:
            data = response.json()
            props = data["PropertyTable"]["Properties"][0]
            return {
                "MolecularWeight": props.get("MolecularWeight", None),
                "XLogP": props.get("XLogP", None)
            }
        else:
            st.error(f"'{compound_name}' 조회 실패 (응답 코드: {response.status_code})")
            return None
    except Exception as e:
        st.error(f"'{compound_name}' 통신 에러: {e}")
        return None

# 데이터 로딩 진행 상태 표시 및 데이터프레임 생성
data_list = []
with st.spinner("PubChem API에서 데이터를 불러오는 중..."):
    for alc in ALCOHOLS:
        api_data = fetch_pubchem_data(alc["name"])
        mw = api_data.get("MolecularWeight") if api_data else None
        xlogp = api_data.get("XLogP") if api_data else None
        
        data_list.append({
            "화학명": alc["name"],
            "한글명": alc["kr_name"],
            "탄소 원자 수": alc["carbons"],
            "분자량 (g/mol)": float(mw) if mw is not None else None,
            "XLogP": float(xlogp) if xlogp is not None else None
        })

df = pd.DataFrame(data_list)

# 1. 물질별 분자량과 XLogP 비교표
st.subheader("1. 알코올 5종 물성 비교표")
st.dataframe(df, use_container_width=True)

# 결측치 체크
valid_df = df.dropna(subset=["분자량 (g/mol)", "XLogP"])

if len(valid_df) < len(df):
    st.warning("⚠️ 일부 물질의 데이터를 가져오지 못해 시각화 및 분석에서 제외될 수 있습니다.")

st.markdown("---")

# 2. 탄소 원자 수와 분자량의 관계 그래프 및 추세선
st.subheader("2. 탄소 원자 수 vs 분자량 (Molecular Weight)")
col1, col2 = st.columns([2, 1])

with col1:
    fig_mw = px.scatter(
        valid_df,
        x="탄소 원자 수",
        y="분자량 (g/mol)",
        text="한글명",
        trendline="ols",
        title="탄소 수 증가에 따른 분자량 변화",
        labels={"탄소 원자 수": "탄소 원자 수 (개)", "분자량 (g/mol)": "분자량 (g/mol)"}
    )
    fig_mw.update_traces(textposition="top center", marker=dict(size=12, color="blue"))
    st.plotly_chart(fig_mw, use_container_width=True)

with col2:
    st.markdown("### 📌 과학적 해석")
    st.write(
        "- 탄소 수가 1개 증가할 때마다 **$\text{CH}_2$ (메틸렌기, 약 14.03 g/mol)** 단위가 추가됩니다.\n"
        "- 이에 따라 탄소 원자 수와 분자량 사이에는 매우 명확한 **선형 비례 관계**가 성립합니다."
    )

st.markdown("---")

# 3. 탄소 원자 수와 XLogP의 관계 그래프 및 상관계수
st.subheader("3. 탄소 원자 수 vs XLogP (소수성 지표)")
col3, col4 = st.columns([2, 1])

# 상관계수 계산
if len(valid_df) >= 2:
    corr, p_value = pearsonr(valid_df["탄소 원자 수"], valid_df["XLogP"])
else:
    corr, p_value = 0, 1

with col3:
    fig_xlogp = px.scatter(
        valid_df,
        x="탄소 원자 수",
        y="XLogP",
        text="한글명",
        trendline="ols",
        title="탄소 수 증가에 따른 XLogP 변화",
        labels={"탄소 원자 수": "탄소 원자 수 (개)", "XLogP": "XLogP"}
    )
    fig_xlogp.update_traces(textposition="top center", marker=dict(size=12, color="orange"))
    st.plotly_chart(fig_xlogp.update_layout(), use_container_width=True)

with col4:
    st.metric(label="피어슨 상관계수 (r)", value=f"{corr:.4f}")
    st.markdown("### 📌 과학적 해석")
    st.write(
        "- **XLogP 정의**: 옥탄올/물 분배계수($\log P$)의 **계산값(Computational Value)**으로, 값이 클수록 **소수성(비극성)**이 강함을 의미합니다.\n"
        "- 탄소 사슬이 길어질수록 극성인 히드록시기(-OH)의 영향에 비해 비극성인 알킬기 사슬의 영향력이 커져 전체적인 **소수성이 증가**합니다."
    )

st.markdown("---")

# 4. 분석 결과 종합 및 유의사항
st.subheader("4. 종합 분석 및 분석 시 유의사항")

st.info(
    "💡 **데이터 출처 안내**: 본 프로그램의 데이터는 **PubChem API**를 통해 실시간 수집되었으며, "
    "실험값을 직접 측정하지 않고 API에서 제공하는 계산치(XLogP)를 활용하였습니다."
)

st.warning(
    "⚠️ **주의사항 (상관관계 vs 인과관계)**\n\n"
    "탄소 원자 수와 XLogP 사이에 높은 양의 상관관계가 나타나지만, '탄소 원자 수가 많다'는 사실 자체가 직접적으로 소수성을 유발하는 유일한 원인은 아닙니다.\n\n"
    "실제 원인은 탄소 수 증가에 따라 **분자 내 비극성 알킬 사슬의 부피가 커지고, 분자 간 런던 분산력이 증가**하기 때문입니다. "
    "따라서 통계적 상관관계를 해석할 때는 물리화학적 메커니즘을 함께 고려해야 합니다."
)

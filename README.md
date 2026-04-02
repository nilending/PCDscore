XGBoost model file for quantifying programmed cell death‑related risk in colorectal cancer.
Load using xgb.load("PCDscore.model") and apply via predict().
Takes transcriptomic feature vector as input and returns continuous PCDscore values.
Compatible with the xgboost R package.

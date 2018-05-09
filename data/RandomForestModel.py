from sklearn.ensemble import RandomForestRegressor

class RandomForestModel:
    """Class for storing random forest model"""
    def __init__(self, regressor: RandomForestRegressor, features_mean: dict,
                 order:list):
        self.regressor = regressor
        self.features_mean = features_mean
        self.order = order

    def predict(self, features):
        return self.regressor.predict(features[self.order])
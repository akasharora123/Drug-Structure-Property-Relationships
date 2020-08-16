
   def best_RidgeCV(self, alphas, cv, score):

       model = RidgeCV(alphas=alphas, cv=cv, scoring=score)
       model = model.fit(self.X, self.y)
       print("Best Ridge model")
       print("Parameters:", model.get_params)
       print("Score: ", -1*model.score(self.X, self.y)) 
       
       return model

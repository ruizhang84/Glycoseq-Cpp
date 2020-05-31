#ifndef ENGINE_ANALYSIS_SVM_ANALYZER_H
#define ENGINE_ANALYSIS_SVM_ANALYZER_H

#include <algorithm>
#include <vector>
#include "../../lib/svm.h"
#include "../search/search_result.h"

namespace engine {
namespace analysis {

using namespace engine::search;

class SVMAnalyzer
{
public:
    SVMAnalyzer()
    {
        model_ = nullptr;
        problem_ = new svm_problem();
        parameter_ = new svm_parameter();

        // set up parameter
        parameter_->svm_type = C_SVC;
        parameter_->kernel_type = POLY;
        parameter_->degree = 3;
        parameter_->C = 1.0; //
        parameter_->gamma = 1.0; //
        parameter_->eps = 0.001;
        parameter_->nu = 0.5;
        parameter_->p = 0.10;
        parameter_->cache_size = 100;
        parameter_->shrinking = 1;
        parameter_->probability = 1;
    }

    ~SVMAnalyzer()
    {
        if (model_ != nullptr)
        {
            ProblemClear();
            ModelClear();	
        }
        else if (problem_->l > 0)
        {
            ProblemClear();
        }
        delete problem_;
        delete parameter_;
    }

    virtual void Training(std::vector<SearchResult> targets, std::vector<SearchResult> decoys)
    {
        // set up problem
        set_problem(targets, decoys);
        
        // training the model
        if (model_ != nullptr)
        {
            ProblemClear();
            ModelClear();	
        }

        model_ = svm_train(problem_, parameter_);          
    }


    virtual double Predicting(const SearchResult& result)
    {
        svm_node* node = SVMNode(result);
        double pred = svm_predict(model_, node);
        delete[] node;
        return pred;
    }

    virtual double PredictingProbability(const SearchResult& result)
    {
        svm_node* node = SVMNode(result);
        double prob = 0;
        svm_predict_probability(model_, node, &prob);
        delete[] node;
        return prob;
    }

    virtual void set_problem(std::vector<SearchResult> targets, std::vector<SearchResult> decoys)
    {
        if (problem_->l > 0)
        {
            ProblemClear();
        }

        int size = (int) targets.size() + decoys.size();
        int index = 0;
        svm_node** nodes_arr = new svm_node*[size];
        double* y_arr = new double[size];

        for(const auto& it : targets)
        {
            *(nodes_arr + index) = SVMNode(it);
            *(y_arr + index) = 1;
            index++;
        }
        for(const auto& it : decoys)
        {
            *(nodes_arr + index) = SVMNode(it);
            *(y_arr + index) = -1;
            index++;
        }
       
        // set up problem
        problem_->l = size;
        problem_->x = nodes_arr;
        problem_->y = y_arr;  
    }

protected:
    virtual void ModelClear()
    {
        int l = model_->l;
        int k = model_->nr_class;
        free(model_->label);
        free(model_->rho);
        if (model_->probA != NULL)
            free(model_->probA);
        if (model_->probB != NULL)
            free(model_->probB);
        free(model_->nSV);
        for(int i = 0; i < l; i++)
        {
            delete (model_->SV)[i];
        }
        free(model_->sv_indices);
        for(int i = 0; i < k; i++)
        {
            free((model_->sv_coef)[i]);
        }
        free(model_->sv_coef);
    }

    virtual void ProblemClear()
    {
        for (int i = 0; i < problem_->l; i++)
        {
            delete[] (problem_->x)[i];
        }
        delete[] problem_->x;
        delete[] problem_->y;
    }

    virtual svm_node* SVMNode(const SearchResult& result)
    {
        svm_node* node = new svm_node[kDims+1];
        for (int i = 0; i < kDims; i++)
        {
            node[i].index = i;
            node[i].value = 0;
        }
        node[kDims].index = -1;

        for (const auto& it : result.Match())
        {
            switch (it.first)
            {
            case SearchType::Core:
                node[0].value = it.second;
                break;
             case SearchType::Branch:
                node[1].value = it.second;
                break;
            case SearchType::Terminal:
                node[2].value = it.second;
                break;
            case SearchType::Peptide:
                node[3].value = it.second;
                break;
            case SearchType::Oxonium:
                node[4].value = it.second;
                break;      
            default:
                break;
            }
        }
        return node;
    }

    static int constexpr kDims = 5; //glycan core, branch, terminal, peptide, oxonium 
    svm_model* model_;
    svm_problem* problem_;
    svm_parameter* parameter_;
};


} // namespace analysis 
} // namespace engine


#endif
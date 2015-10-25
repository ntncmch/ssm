{% extends "ordered.tpl" %}

{% block code %}

void ssm_compute_diff(ssm_X_t *p_X, double t, ssm_par_t *par, ssm_nav_t *nav, ssm_calc_t *calc)
{
    {% if diff.n_diffs %}

    int i;
    int n_browns = {{ diff.n_browns}};
    ssm_it_states_t *it = nav->states_diff;

    double diffed[it->length];
    double *X = p_X->proj;

    for(i=0; i<it->length; i++){
        ssm_state_t *p = it->p[i];
        diffed[i] = p->f_inv(X[p->offset]);
    }

    double dt = p_X->dt;
    double sqrt_dt = sqrt(dt);
    
    double _w[n_browns];
    for(i=0; i<n_browns; i++){
        _w[i] = gsl_ran_ugaussian(calc->randgsl);
    }

    {% for i in range(diff.n_diffs) %}
    p_X->proj[it->p[{{ loop.index0 }}]->offset] += dt*({{ diff.drift_terms[i] }}) + sqrt_dt*({{ diff.dispersion_terms[i] }});
    {% endfor %}

    {% endif %}
}


{% endblock %}


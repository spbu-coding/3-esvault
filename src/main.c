#include <stdio.h>
#include <math.h>

#define eps_float 1e-6f
#define eps_double 1e-6

typedef struct {
    float x;
    float y;
} float_point_t;

typedef struct {
    double x;
    double y;
} double_point_t;

float_point_t resolve_float_system(float delta);
double_point_t resolve_double_system(double delta);

float calculate_float_dist(float_point_t p1, float_point_t p2);
double calculate_double_dist(double_point_t p1, double_point_t p2);


int main() {
    float delta_float = 0.00001f;
    float_point_t float_result = resolve_float_system(0.0f);
    float_point_t float_result_delta = resolve_float_system(delta_float);
    float float_res_deviation = calculate_float_dist(float_result, float_result_delta);
    printf("==========Results for floats===========\n");
    while (float_res_deviation >= eps_float) {
        printf("Delta: %.8f, deviation: %.8f\n", delta_float, float_res_deviation);
        delta_float /= 2.0f;
        float_result_delta = resolve_float_system(delta_float);
        float_res_deviation = calculate_float_dist(float_result, float_result_delta);
    }
    double delta_double = 0.00001;
    double_point_t double_result = resolve_double_system(0.0);
    double_point_t double_result_delta = resolve_double_system(delta_double);
    double double_res_deviation = calculate_double_dist(double_result, double_result_delta);
    printf("==========Results for doubles===========\n");
    while (double_res_deviation >= eps_double) {
        printf("Delta: %.11f, deviation: %f\n", delta_double, double_res_deviation);
        delta_double /= 2.0;
        double_result_delta = resolve_double_system(delta_double);
        double_res_deviation = calculate_double_dist(double_result, double_result_delta);
    }
    return 0;
}

float_point_t resolve_float_system(float delta) {
    float_point_t result;
    result.y = ((2.00001f + delta) - 2.0f) / 0.0001f;
    result.x = 2.0f - result.y;
    return result;
}

double_point_t resolve_double_system(double delta) {
    double_point_t result;
    result.y = ((2.00001 + delta) - 2.0) / 0.0001;
    result.x = 2.0 - result.y;
    return result;
}

float calculate_float_dist(float_point_t p1, float_point_t p2) {
    return sqrtf(powf(p1.x - p2.x, 2) + powf(p1.y - p2.y, 2));
}

double calculate_double_dist(double_point_t p1, double_point_t p2) {
    return (sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2)));
}

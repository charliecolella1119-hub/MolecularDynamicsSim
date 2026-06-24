#include <SFML/Graphics.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <vector>

// --------------------
// Window / Layout
// --------------------

const int WINDOW_WIDTH = 1000;
const int WINDOW_HEIGHT = 800;

const float BOX_LEFT = 260.0f;
const float BOX_TOP = 70.0f;
const float BOX_WIDTH = 700.0f;
const float BOX_HEIGHT = 650.0f;

const float PARTICLE_RADIUS = 4.0f;
const float SIDEBAR_GRAPH_LABEL_TOP = 605.0f;
const float SIDEBAR_GRAPH_TOP = 630.0f;
const float SIDEBAR_GRAPH_HEIGHT = 95.0f;
const float PI = 3.14159265359f;

// --------------------
// Physics Constants
// --------------------

const float EPSILON = 1.0f;
const float SIGMA = 12.0f;
const float MASS = 1.0f;

const float TIME_STEP = 0.0025f;
const float FORCE_SCALE = 120.0f;

enum class Preset
{
    Gas,
    Liquid,
    Crystal,
    HotGas,
    DenseFluid,
    MeltingCrystal
};

enum class GraphMode
{
    Temperature,
    RadialDistribution,
    VelocityDistribution,
    PhaseDiagram
};

struct Particle
{
    sf::Vector2f position;
    sf::Vector2f velocity;
    sf::Vector2f force;
    std::deque<sf::Vector2f> trail;
};

struct SimulationSample
{
    float time;
    float temperature;
    float pressure;
    float kineticEnergy;
    float potentialEnergy;
    float totalEnergy;
};

float lengthSquared(const sf::Vector2f& v)
{
    return v.x * v.x + v.y * v.y;
}

float speed(const Particle& p)
{
    return std::sqrt(lengthSquared(p.velocity));
}

std::string formatFloat(float value, int precision = 2)
{
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(precision) << value;
    return stream.str();
}

sf::Color speedColor(float v)
{
    if (v < 40.0f)
        return sf::Color(80, 180, 255);

    if (v < 90.0f)
        return sf::Color(80, 255, 140);

    if (v < 150.0f)
        return sf::Color(255, 180, 60);

    return sf::Color(255, 80, 60);
}

void applyPeriodicBoundaries(Particle& p)
{
    if (p.position.x < BOX_LEFT)
        p.position.x += BOX_WIDTH;

    if (p.position.x > BOX_LEFT + BOX_WIDTH)
        p.position.x -= BOX_WIDTH;

    if (p.position.y < BOX_TOP)
        p.position.y += BOX_HEIGHT;

    if (p.position.y > BOX_TOP + BOX_HEIGHT)
        p.position.y -= BOX_HEIGHT;
}

sf::Vector2f minimumImage(sf::Vector2f r)
{
    if (r.x > 0.5f * BOX_WIDTH)
        r.x -= BOX_WIDTH;

    if (r.x < -0.5f * BOX_WIDTH)
        r.x += BOX_WIDTH;

    if (r.y > 0.5f * BOX_HEIGHT)
        r.y -= BOX_HEIGHT;

    if (r.y < -0.5f * BOX_HEIGHT)
        r.y += BOX_HEIGHT;

    return r;
}

void computeForces(std::vector<Particle>& particles)
{
    for (auto& p : particles)
    {
        p.force = sf::Vector2f(0.0f, 0.0f);
    }

    for (size_t i = 0; i < particles.size(); i++)
    {
        for (size_t j = i + 1; j < particles.size(); j++)
        {
            sf::Vector2f r =
                minimumImage(particles[j].position - particles[i].position);

            float r2 = lengthSquared(r);

            if (r2 < 1.0f)
                r2 = 1.0f;

            float invR2 = 1.0f / r2;
            float sigma2OverR2 = (SIGMA * SIGMA) * invR2;

            float sr6 = sigma2OverR2 * sigma2OverR2 * sigma2OverR2;
            float sr12 = sr6 * sr6;

            float forceMagnitude =
                24.0f * EPSILON * (2.0f * sr12 - sr6) * invR2;

            sf::Vector2f force = FORCE_SCALE * forceMagnitude * r;

            particles[i].force -= force;
            particles[j].force += force;
        }
    }
}

void velocityVerletStep(std::vector<Particle>& particles)
{
    for (auto& p : particles)
    {
        p.velocity += 0.5f * TIME_STEP * p.force / MASS;
        p.position += TIME_STEP * p.velocity;
        applyPeriodicBoundaries(p);
    }

    computeForces(particles);

    for (auto& p : particles)
    {
        p.velocity += 0.5f * TIME_STEP * p.force / MASS;
    }
}

void scaleVelocities(std::vector<Particle>& particles, float factor)
{
    for (auto& p : particles)
    {
        p.velocity *= factor;
    }
}

void updateParticleTrails(std::vector<Particle>& particles,
                          size_t maxTrailPoints)
{
    for (auto& p : particles)
    {
        if (!p.trail.empty() &&
            lengthSquared(p.position - p.trail.back()) > 0.20f * BOX_WIDTH * BOX_WIDTH)
        {
            p.trail.clear();
        }

        p.trail.push_back(p.position);

        while (p.trail.size() > maxTrailPoints)
        {
            p.trail.pop_front();
        }
    }
}

float computeKineticEnergy(const std::vector<Particle>& particles)
{
    float kineticEnergy = 0.0f;

    for (const auto& p : particles)
    {
        kineticEnergy += 0.5f * MASS * lengthSquared(p.velocity);
    }

    return kineticEnergy;
}

float computePressure(const std::vector<Particle>& particles,
                      float temperature)
{
    float virial = 0.0f;

    for (size_t i = 0; i < particles.size(); i++)
    {
        for (size_t j = i + 1; j < particles.size(); j++)
        {
            sf::Vector2f r =
                minimumImage(particles[j].position - particles[i].position);

            float r2 = lengthSquared(r);

            if (r2 < 1.0f)
                r2 = 1.0f;

            float invR2 = 1.0f / r2;
            float sigma2OverR2 = (SIGMA * SIGMA) * invR2;

            float sr6 = sigma2OverR2 * sigma2OverR2 * sigma2OverR2;
            float sr12 = sr6 * sr6;

            float forceMagnitude =
                24.0f * EPSILON * (2.0f * sr12 - sr6) * invR2;

            sf::Vector2f force = FORCE_SCALE * forceMagnitude * r;
            virial += r.x * force.x + r.y * force.y;
        }
    }

    float area = BOX_WIDTH * BOX_HEIGHT;
    float idealPressure =
        static_cast<float>(particles.size()) * temperature / area;

    return idealPressure + 0.5f * virial / area;
}

float computePotentialEnergy(const std::vector<Particle>& particles)
{
    float potentialEnergy = 0.0f;

    for (size_t i = 0; i < particles.size(); i++)
    {
        for (size_t j = i + 1; j < particles.size(); j++)
        {
            sf::Vector2f r =
                minimumImage(particles[j].position - particles[i].position);

            float r2 = lengthSquared(r);

            if (r2 < 1.0f)
                r2 = 1.0f;

            float invR2 = 1.0f / r2;
            float sigma2OverR2 = (SIGMA * SIGMA) * invR2;

            float sr6 = sigma2OverR2 * sigma2OverR2 * sigma2OverR2;
            float sr12 = sr6 * sr6;

            potentialEnergy += 4.0f * EPSILON * (sr12 - sr6);
        }
    }

    return FORCE_SCALE * potentialEnergy;
}

int particleCountForPreset(Preset preset)
{
    if (preset == Preset::DenseFluid)
        return 196;

    return 144;
}

std::vector<Particle> createParticles(int count, Preset preset)
{
    std::vector<Particle> particles;

    std::mt19937 rng(static_cast<unsigned int>(std::random_device{}()));

    float velocityRange = 35.0f;

    if (preset == Preset::Gas)
        velocityRange = 45.0f;
    else if (preset == Preset::Liquid)
        velocityRange = 25.0f;
    else if (preset == Preset::Crystal)
        velocityRange = 6.0f;
    else if (preset == Preset::HotGas)
        velocityRange = 85.0f;
    else if (preset == Preset::DenseFluid)
        velocityRange = 18.0f;
    else if (preset == Preset::MeltingCrystal)
        velocityRange = 28.0f;

    std::uniform_real_distribution<float> velocityDist(
        -velocityRange,
        velocityRange
    );

    int cols = static_cast<int>(std::sqrt(count));
    int rows = cols;

    float spacingX = BOX_WIDTH / static_cast<float>(cols + 1);
    float spacingY = BOX_HEIGHT / static_cast<float>(rows + 1);

    for (int y = 0; y < rows; y++)
    {
        for (int x = 0; x < cols; x++)
        {
            if (static_cast<int>(particles.size()) >= count)
                break;

            Particle p;

            if (preset == Preset::Gas || preset == Preset::HotGas)
            {
                p.position = sf::Vector2f(
                    BOX_LEFT + spacingX * (x + 1),
                    BOX_TOP + spacingY * (y + 1)
                );
            }
            else if (preset == Preset::Liquid || preset == Preset::DenseFluid)
            {
                float clusterWidth = 320.0f;
                float clusterHeight = 260.0f;

                if (preset == Preset::DenseFluid)
                {
                    clusterWidth = 420.0f;
                    clusterHeight = 360.0f;
                }

                p.position = sf::Vector2f(
                    BOX_LEFT + BOX_WIDTH * 0.5f -
                        clusterWidth * 0.5f +
                        (x + 1) *
                        (clusterWidth / static_cast<float>(cols + 1)),

                    BOX_TOP + BOX_HEIGHT * 0.5f -
                        clusterHeight * 0.5f +
                        (y + 1) *
                        (clusterHeight / static_cast<float>(rows + 1))
                );
            }
            else
            {
                float crystalSpacing = 22.0f;

                if (preset == Preset::MeltingCrystal)
                    crystalSpacing = 20.0f;

                p.position = sf::Vector2f(
                    BOX_LEFT + BOX_WIDTH * 0.5f -
                        cols * crystalSpacing * 0.5f +
                        x * crystalSpacing,

                    BOX_TOP + BOX_HEIGHT * 0.5f -
                        rows * crystalSpacing * 0.5f +
                        y * crystalSpacing
                );
            }

            p.velocity = sf::Vector2f(
                velocityDist(rng),
                velocityDist(rng)
            );

            p.force = sf::Vector2f(0.0f, 0.0f);
            p.trail.push_back(p.position);

            particles.push_back(p);
        }
    }

    return particles;
}

void updateTemperatureHistory(std::vector<float>& history,
                              float temperature,
                              size_t maxPoints)
{
    history.push_back(temperature);

    if (history.size() > maxPoints)
    {
        history.erase(history.begin());
    }
}

void drawTemperatureGraph(sf::RenderWindow& window,
                          const std::vector<float>& history,
                          sf::Vector2f position,
                          sf::Vector2f size)
{
    sf::RectangleShape graphBox;
    graphBox.setPosition(position);
    graphBox.setSize(size);
    graphBox.setFillColor(sf::Color(18, 21, 30));
    graphBox.setOutlineColor(sf::Color(70, 80, 100));
    graphBox.setOutlineThickness(1.0f);

    window.draw(graphBox);

    if (history.size() < 2)
        return;

    float minTemp = history[0];
    float maxTemp = history[0];

    for (float t : history)
    {
        if (t < minTemp) minTemp = t;
        if (t > maxTemp) maxTemp = t;
    }

    if (std::abs(maxTemp - minTemp) < 0.001f)
    {
        maxTemp = minTemp + 1.0f;
    }

    sf::VertexArray line(sf::PrimitiveType::LineStrip);

    for (size_t i = 0; i < history.size(); i++)
    {
        float x =
            position.x +
            (static_cast<float>(i) /
             static_cast<float>(history.size() - 1)) *
            size.x;

        float normalized =
            (history[i] - minTemp) / (maxTemp - minTemp);

        float y =
            position.y +
            size.y -
            normalized * size.y;

        sf::Vertex vertex;
        vertex.position = sf::Vector2f(x, y);
        vertex.color = sf::Color(80, 220, 255);

        line.append(vertex);
    }

    window.draw(line);
}

std::string presetName(Preset preset)
{
    if (preset == Preset::Gas)
        return "Gas";

    if (preset == Preset::Liquid)
        return "Liquid";

    if (preset == Preset::Crystal)
        return "Crystal";

    if (preset == Preset::HotGas)
        return "Hot Gas";

    if (preset == Preset::DenseFluid)
        return "Dense Fluid";

    return "Melting Crystal";
}

bool saveThermodynamicData(const std::vector<SimulationSample>& samples,
                           const std::string& path)
{
    std::ofstream file(path);

    if (!file)
        return false;

    file << "time,temperature,pressure,kinetic_energy,potential_energy,total_energy\n";

    for (const auto& sample : samples)
    {
        file << sample.time << ','
             << sample.temperature << ','
             << sample.pressure << ','
             << sample.kineticEnergy << ','
             << sample.potentialEnergy << ','
             << sample.totalEnergy << '\n';
    }

    return true;
}

std::vector<float> computeRadialDistribution(
    const std::vector<Particle>& particles,
    int binCount,
    float maxRadius)
{
    std::vector<float> bins(binCount, 0.0f);

    float binWidth = maxRadius / static_cast<float>(binCount);

    for (size_t i = 0; i < particles.size(); i++)
    {
        for (size_t j = i + 1; j < particles.size(); j++)
        {
            sf::Vector2f r =
                minimumImage(particles[j].position - particles[i].position);

            float distance = std::sqrt(lengthSquared(r));

            if (distance < maxRadius)
            {
                int bin = static_cast<int>(distance / binWidth);

                if (bin >= 0 && bin < binCount)
                {
                    bins[bin] += 2.0f;
                }
            }
        }
    }

    float density =
        static_cast<float>(particles.size()) / (BOX_WIDTH * BOX_HEIGHT);

    for (int i = 0; i < binCount; i++)
    {
        float r1 = i * binWidth;
        float r2 = r1 + binWidth;

        float shellArea = PI * (r2 * r2 - r1 * r1);

        float idealCount =
            density * shellArea *
            static_cast<float>(particles.size());

        if (idealCount > 0.0f)
        {
            bins[i] /= idealCount;
        }
    }

    return bins;
}

void drawBarGraph(sf::RenderWindow& window,
                  const std::vector<float>& values,
                  sf::Vector2f position,
                  sf::Vector2f size,
                  sf::Color barColor)
{
    sf::RectangleShape graphBox;
    graphBox.setPosition(position);
    graphBox.setSize(size);
    graphBox.setFillColor(sf::Color(18, 21, 30));
    graphBox.setOutlineColor(sf::Color(70, 80, 100));
    graphBox.setOutlineThickness(1.0f);

    window.draw(graphBox);

    if (values.empty())
        return;

    float maxValue = 0.0f;

    for (float v : values)
    {
        if (v > maxValue)
            maxValue = v;
    }

    if (maxValue < 0.001f)
        maxValue = 1.0f;

    float barWidth = size.x / static_cast<float>(values.size());

    for (size_t i = 0; i < values.size(); i++)
    {
        float normalized = values[i] / maxValue;

        if (normalized > 1.0f)
            normalized = 1.0f;

        float barHeight = normalized * size.y;

        sf::RectangleShape bar;
        bar.setPosition(sf::Vector2f(
            position.x + i * barWidth,
            position.y + size.y - barHeight
        ));

        bar.setSize(sf::Vector2f(
            barWidth - 1.0f,
            barHeight
        ));

        bar.setFillColor(barColor);

        window.draw(bar);
    }
}

std::vector<float> computeVelocityHistogram(
    const std::vector<Particle>& particles,
    int bins,
    float maxVelocity)
{
    std::vector<float> histogram(bins, 0.0f);

    for (const auto& p : particles)
    {
        float speed = std::sqrt(p.velocity.x * p.velocity.x +
                                p.velocity.y * p.velocity.y);

        int bin =
            static_cast<int>((speed / maxVelocity) * bins);

        if (bin >= 0 && bin < bins)
            histogram[bin] += 1.0f;
    }

    return histogram;
}


void drawRadialDistributionGraph(sf::RenderWindow& window,
                                 const std::vector<float>& values,
                                 sf::Vector2f position,
                                 sf::Vector2f size)
{
    sf::RectangleShape graphBox;
    graphBox.setPosition(position);
    graphBox.setSize(size);
    graphBox.setFillColor(sf::Color(18, 21, 30));
    graphBox.setOutlineColor(sf::Color(70, 80, 100));
    graphBox.setOutlineThickness(1.0f);
    window.draw(graphBox);

    if (values.empty())
        return;

    float maxValue = std::max(2.0f, *std::max_element(values.begin(), values.end()));
    float baselineY = position.y + size.y - (1.0f / maxValue) * size.y;

    sf::VertexArray idealLine(sf::PrimitiveType::Lines);
    sf::Vertex idealStart;
    idealStart.position = sf::Vector2f(position.x, baselineY);
    idealStart.color = sf::Color(90, 105, 125);
    sf::Vertex idealEnd;
    idealEnd.position = sf::Vector2f(position.x + size.x, baselineY);
    idealEnd.color = sf::Color(90, 105, 125);
    idealLine.append(idealStart);
    idealLine.append(idealEnd);
    window.draw(idealLine);

    float barWidth = size.x / static_cast<float>(values.size());

    for (size_t i = 0; i < values.size(); i++)
    {
        float normalized = std::clamp(values[i] / maxValue, 0.0f, 1.0f);
        float barHeight = normalized * size.y;

        sf::RectangleShape bar;
        bar.setPosition(sf::Vector2f(
            position.x + i * barWidth,
            position.y + size.y - barHeight
        ));
        bar.setSize(sf::Vector2f(std::max(1.0f, barWidth - 1.0f), barHeight));
        bar.setFillColor(sf::Color(255, 180, 60));
        window.draw(bar);
    }
}

void drawVelocityDistributionGraph(sf::RenderWindow& window,
                                   const std::vector<float>& histogram,
                                   float temperature,
                                   float maxVelocity,
                                   sf::Vector2f position,
                                   sf::Vector2f size)
{
    sf::RectangleShape graphBox;
    graphBox.setPosition(position);
    graphBox.setSize(size);
    graphBox.setFillColor(sf::Color(18, 21, 30));
    graphBox.setOutlineColor(sf::Color(70, 80, 100));
    graphBox.setOutlineThickness(1.0f);
    window.draw(graphBox);

    if (histogram.empty())
        return;

    float totalCount = 0.0f;
    float maxHistogram = 0.0f;

    for (float value : histogram)
    {
        totalCount += value;
        maxHistogram = std::max(maxHistogram, value);
    }

    float binWidth = maxVelocity / static_cast<float>(histogram.size());
    std::vector<float> fit(histogram.size(), 0.0f);
    float maxFit = 0.0f;

    if (temperature > 0.001f)
    {
        for (size_t i = 0; i < fit.size(); i++)
        {
            float v = (static_cast<float>(i) + 0.5f) * binWidth;
            float density =
                (MASS / temperature) * v *
                std::exp(-(MASS * v * v) / (2.0f * temperature));
            fit[i] = totalCount * density * binWidth;
            maxFit = std::max(maxFit, fit[i]);
        }
    }

    float graphMax = std::max({1.0f, maxHistogram, maxFit});
    float barWidth = size.x / static_cast<float>(histogram.size());

    for (size_t i = 0; i < histogram.size(); i++)
    {
        float barHeight = (histogram[i] / graphMax) * size.y;

        sf::RectangleShape bar;
        bar.setPosition(sf::Vector2f(
            position.x + i * barWidth,
            position.y + size.y - barHeight
        ));
        bar.setSize(sf::Vector2f(std::max(1.0f, barWidth - 1.0f), barHeight));
        bar.setFillColor(sf::Color(80, 220, 255));
        window.draw(bar);
    }

    sf::VertexArray fitLine(sf::PrimitiveType::LineStrip);

    for (size_t i = 0; i < fit.size(); i++)
    {
        float x =
            position.x +
            (static_cast<float>(i) + 0.5f) *
            size.x / static_cast<float>(fit.size());
        float y = position.y + size.y - (fit[i] / graphMax) * size.y;

        sf::Vertex vertex;
        vertex.position = sf::Vector2f(x, y);
        vertex.color = sf::Color(255, 95, 170);
        fitLine.append(vertex);
    }

    window.draw(fitLine);
}

void drawPhaseDiagram(sf::RenderWindow& window,
                      const std::vector<sf::Vector2f>& samples,
                      sf::Vector2f currentPoint,
                      sf::Vector2f position,
                      sf::Vector2f size)
{
    sf::RectangleShape graphBox;
    graphBox.setPosition(position);
    graphBox.setSize(size);
    graphBox.setFillColor(sf::Color(18, 21, 30));
    graphBox.setOutlineColor(sf::Color(70, 80, 100));
    graphBox.setOutlineThickness(1.0f);
    window.draw(graphBox);

    if (samples.empty())
        return;

    float minTemp = currentPoint.x;
    float maxTemp = currentPoint.x;
    float minPressure = currentPoint.y;
    float maxPressure = currentPoint.y;

    for (const auto& sample : samples)
    {
        minTemp = std::min(minTemp, sample.x);
        maxTemp = std::max(maxTemp, sample.x);
        minPressure = std::min(minPressure, sample.y);
        maxPressure = std::max(maxPressure, sample.y);
    }

    if (std::abs(maxTemp - minTemp) < 0.001f)
        maxTemp = minTemp + 1.0f;

    if (std::abs(maxPressure - minPressure) < 0.001f)
        maxPressure = minPressure + 1.0f;

    sf::CircleShape point(2.0f);
    point.setOrigin(sf::Vector2f(2.0f, 2.0f));

    for (size_t i = 0; i < samples.size(); i++)
    {
        float xNorm = (samples[i].x - minTemp) / (maxTemp - minTemp);
        float yNorm = (samples[i].y - minPressure) / (maxPressure - minPressure);

        point.setPosition(sf::Vector2f(
            position.x + xNorm * size.x,
            position.y + size.y - yNorm * size.y
        ));

        std::uint8_t alpha =
            static_cast<std::uint8_t>(60 + 180 * (static_cast<float>(i + 1) /
                                               static_cast<float>(samples.size())));
        point.setFillColor(sf::Color(115, 230, 155, alpha));
        window.draw(point);
    }

    point.setRadius(3.5f);
    point.setOrigin(sf::Vector2f(3.5f, 3.5f));
    point.setFillColor(sf::Color(255, 235, 120));

    float xNorm = (currentPoint.x - minTemp) / (maxTemp - minTemp);
    float yNorm = (currentPoint.y - minPressure) / (maxPressure - minPressure);

    point.setPosition(sf::Vector2f(
        position.x + xNorm * size.x,
        position.y + size.y - yNorm * size.y
    ));
    window.draw(point);
}

void drawParticleTrails(sf::RenderWindow& window,
                        const std::vector<Particle>& particles)
{
    for (const auto& p : particles)
    {
        if (p.trail.size() < 2)
            continue;

        sf::VertexArray trail(sf::PrimitiveType::LineStrip);

        for (size_t i = 0; i < p.trail.size(); i++)
        {
            if (i > 0 &&
                lengthSquared(p.trail[i] - p.trail[i - 1]) >
                    0.15f * BOX_WIDTH * BOX_WIDTH)
            {
                if (trail.getVertexCount() > 1)
                    window.draw(trail);

                trail.clear();
            }

            sf::Vertex vertex;
            vertex.position = p.trail[i];
            std::uint8_t alpha =
                static_cast<std::uint8_t>(25 + 115 *
                    (static_cast<float>(i + 1) /
                     static_cast<float>(p.trail.size())));
            vertex.color = sf::Color(90, 180, 255, alpha);
            trail.append(vertex);
        }

        if (trail.getVertexCount() > 1)
            window.draw(trail);
    }
}

int main()
{
    sf::RenderWindow window(
        sf::VideoMode(sf::Vector2u(WINDOW_WIDTH, WINDOW_HEIGHT)),
        "Molecular Dynamics Simulator"
    );

    window.setFramerateLimit(60);

    sf::Font font;

    if (!font.openFromFile("../assets/Arial.ttf"))
    {
        std::cerr << "Failed to load font\n";
    }

    sf::Text infoText(font);
    infoText.setCharacterSize(15);
    infoText.setFillColor(sf::Color(230, 230, 230));
    infoText.setPosition(sf::Vector2f(20.0f, 85.0f));

    sf::Text graphLabel(font);
    graphLabel.setCharacterSize(14);
    graphLabel.setFillColor(sf::Color(200, 220, 230));
    graphLabel.setPosition(sf::Vector2f(25.0f, SIDEBAR_GRAPH_LABEL_TOP));
    graphLabel.setString("Temperature History");

    sf::RectangleShape uiPanel;
    uiPanel.setPosition(sf::Vector2f(10.0f, 70.0f));
    uiPanel.setSize(sf::Vector2f(230.0f, 650.0f));
    uiPanel.setFillColor(sf::Color(25, 28, 38));
    uiPanel.setOutlineColor(sf::Color(70, 80, 100));
    uiPanel.setOutlineThickness(1.0f);

    sf::RectangleShape box;
    box.setPosition(sf::Vector2f(BOX_LEFT, BOX_TOP));
    box.setSize(sf::Vector2f(BOX_WIDTH, BOX_HEIGHT));
    box.setFillColor(sf::Color::Transparent);
    box.setOutlineColor(sf::Color(180, 180, 180));
    box.setOutlineThickness(2.0f);

    sf::CircleShape particleShape(PARTICLE_RADIUS);
    particleShape.setOrigin(sf::Vector2f(PARTICLE_RADIUS, PARTICLE_RADIUS));

    Preset currentPreset = Preset::Gas;
    std::vector<Particle> particles = createParticles(
        particleCountForPreset(currentPreset),
        currentPreset
    );
    computeForces(particles);

    bool paused = false;
    bool showTrails = true;
    bool csvSaved = false;
    float simulationTime = 0.0f;

    std::vector<float> temperatureHistory;
    const size_t maxHistoryPoints = 180;

    float smoothedTemperature = 0.0f;
    bool hasSmoothedTemperature = false;

    GraphMode graphMode = GraphMode::Temperature;

    std::vector<float> radialDistribution;
    const int radialBinCount = 45;
    const float radialMaxDistance = 180.0f;

    int radialFrameCounter = 0;

    std::vector<float> velocityHistogram;

    const int velocityBinCount = 30;
    float velocityMax = 300.0f;

    std::vector<sf::Vector2f> phaseSamples;
    const size_t maxPhaseSamples = 240;
    int phaseFrameCounter = 0;
    const size_t maxTrailPoints = 28;

    std::vector<SimulationSample> thermodynamicLog;
    int logFrameCounter = 0;
    const size_t maxLoggedSamples = 2000;

    auto resetSimulation = [&](Preset preset)
    {
        currentPreset = preset;
        particles = createParticles(particleCountForPreset(currentPreset), currentPreset);
        computeForces(particles);
        temperatureHistory.clear();
        phaseSamples.clear();
        radialDistribution.clear();
        velocityHistogram.clear();
        thermodynamicLog.clear();
        radialFrameCounter = 0;
        phaseFrameCounter = 0;
        logFrameCounter = 0;
        simulationTime = 0.0f;
        csvSaved = false;
        hasSmoothedTemperature = false;
    };

    while (window.isOpen())
    {
        while (const std::optional event = window.pollEvent())
        {
            if (event->is<sf::Event::Closed>())
            {
                window.close();
            }

            if (const auto* keyPressed =
                    event->getIf<sf::Event::KeyPressed>())
            {
                if (keyPressed->code == sf::Keyboard::Key::Escape)
                {
                    window.close();
                }

                if (keyPressed->code == sf::Keyboard::Key::Space)
                {
                    paused = !paused;
                }

                if (keyPressed->code == sf::Keyboard::Key::R)
                {
                    resetSimulation(currentPreset);
                }

                if (keyPressed->code == sf::Keyboard::Key::Up)
                {
                    scaleVelocities(particles, 1.10f);
                }

                if (keyPressed->code == sf::Keyboard::Key::Down)
                {
                    scaleVelocities(particles, 0.90f);
                }

                if (keyPressed->code == sf::Keyboard::Key::Num1)
                {
                    resetSimulation(Preset::Gas);
                }

                if (keyPressed->code == sf::Keyboard::Key::Num2)
                {
                    resetSimulation(Preset::Liquid);
                }

                if (keyPressed->code == sf::Keyboard::Key::Num3)
                {
                    resetSimulation(Preset::Crystal);
                }

                if (keyPressed->code == sf::Keyboard::Key::Num4)
                {
                    resetSimulation(Preset::HotGas);
                }

                if (keyPressed->code == sf::Keyboard::Key::Num5)
                {
                    resetSimulation(Preset::DenseFluid);
                }

                if (keyPressed->code == sf::Keyboard::Key::Num6)
                {
                    resetSimulation(Preset::MeltingCrystal);
                }

                if (keyPressed->code == sf::Keyboard::Key::H)
                {
                    showTrails = !showTrails;
                }

                if (keyPressed->code == sf::Keyboard::Key::C)
                {
                    csvSaved = saveThermodynamicData(
                        thermodynamicLog,
                        "thermodynamic_log.csv"
                    );
                }

                if (keyPressed->code == sf::Keyboard::Key::T)
                {
                    graphMode = GraphMode::Temperature;
                }

                if (keyPressed->code == sf::Keyboard::Key::G)
                {
                    graphMode = GraphMode::RadialDistribution;
                }

                if (keyPressed->code == sf::Keyboard::Key::V)
                {
                    graphMode = GraphMode::VelocityDistribution;
                }

                if (keyPressed->code == sf::Keyboard::Key::P)
                {
                    graphMode = GraphMode::PhaseDiagram;
                }
            }
        }

        if (!paused)
        {
            for (int i = 0; i < 4; i++)
            {
                velocityVerletStep(particles);
                simulationTime += TIME_STEP;
            }

            if (showTrails)
            {
                updateParticleTrails(particles, maxTrailPoints);
            }
        }

        float kineticEnergy = computeKineticEnergy(particles);
        float potentialEnergy = computePotentialEnergy(particles);
        float totalEnergy = kineticEnergy + potentialEnergy;
        float temperature =
            kineticEnergy / static_cast<float>(particles.size());
        float pressure = computePressure(particles, temperature);

        if (!paused)
        {
            if (!hasSmoothedTemperature)
            {
                smoothedTemperature = temperature;
                hasSmoothedTemperature = true;
            }

            smoothedTemperature =
                0.95f * smoothedTemperature +
                0.05f * temperature;

            updateTemperatureHistory(
                temperatureHistory,
                smoothedTemperature,
                maxHistoryPoints
            );
        }

        if (!paused)
        {
            radialFrameCounter++;

            if (radialFrameCounter % 15 == 0)
            {
                radialDistribution = computeRadialDistribution(
                    particles,
                    radialBinCount,
                    radialMaxDistance
                );
            }
        }
        
        velocityHistogram =
            computeVelocityHistogram(
                particles,
                velocityBinCount,
                velocityMax
            );

        if (!paused)
        {
            phaseFrameCounter++;

            if (phaseFrameCounter % 10 == 0)
            {
                phaseSamples.push_back(sf::Vector2f(temperature, pressure));

                if (phaseSamples.size() > maxPhaseSamples)
                {
                    phaseSamples.erase(phaseSamples.begin());
                }
            }

            logFrameCounter++;

            if (logFrameCounter % 10 == 0)
            {
                thermodynamicLog.push_back(SimulationSample{
                    simulationTime,
                    temperature,
                    pressure,
                    kineticEnergy,
                    potentialEnergy,
                    totalEnergy
                });

                if (thermodynamicLog.size() > maxLoggedSamples)
                {
                    thermodynamicLog.erase(thermodynamicLog.begin());
                }
            }
        }

        infoText.setString(
            "Molecular Dynamics Simulator\n"
            "Preset: " + presetName(currentPreset) + "\n"
            "Temp: " + formatFloat(temperature) + "\n" +
            "Pressure: " + formatFloat(pressure, 4) + "\n" +
            "Kinetic Energy: " + formatFloat(kineticEnergy) + "\n" +
            "Potential Energy: " + formatFloat(potentialEnergy) + "\n" +
            "Total Energy: " + formatFloat(totalEnergy) + "\n\n" +
            "Controls:\n"
            "1 Gas  2 Liquid  3 Crystal\n"
            "4 Hot  5 Dense  6 Melt\n"
            "Up/Down = Heat/Cool\n"
            "T/G/V/P = Graph Modes\n"
            "H = Trails " + std::string(showTrails ? "On" : "Off") + "\n" +
            "C = Export CSV" + std::string(csvSaved ? " Saved" : "") + "\n" +
            "Space = Pause\n"
            "R = Reset\n"
            "Esc = Quit"
        );

        window.setTitle(
            "MD Sim | Preset: " + presetName(currentPreset) +
            " | Temp: " + formatFloat(temperature) +
            " | Pressure: " + formatFloat(pressure, 4)
        );

        window.clear(sf::Color(15, 17, 23));

        window.draw(box);

        if (showTrails)
        {
            drawParticleTrails(window, particles);
        }

        for (const auto& p : particles)
        {
            particleShape.setPosition(p.position);
            particleShape.setFillColor(speedColor(speed(p)));
            window.draw(particleShape);
        }

        window.draw(uiPanel);
        window.draw(infoText);

        if (graphMode == GraphMode::Temperature)
        {
            graphLabel.setString("Temperature History");
            window.draw(graphLabel);

            drawTemperatureGraph(
                window,
                temperatureHistory,
                sf::Vector2f(20.0f, SIDEBAR_GRAPH_TOP),
                sf::Vector2f(210.0f, SIDEBAR_GRAPH_HEIGHT)
            );
        }
        else if (graphMode == GraphMode::RadialDistribution)
        {
            graphLabel.setString("Radial Distribution g(r)");
            window.draw(graphLabel);

            drawRadialDistributionGraph(
                window,
                radialDistribution,
                sf::Vector2f(20.0f, SIDEBAR_GRAPH_TOP),
                sf::Vector2f(210.0f, SIDEBAR_GRAPH_HEIGHT)
            );
        }
        else if (graphMode == GraphMode::VelocityDistribution)
        {
            graphLabel.setString("Velocity Histogram");
            window.draw(graphLabel);

            drawVelocityDistributionGraph(
                window,
                velocityHistogram,
                temperature,
                velocityMax,
                sf::Vector2f(20.0f, SIDEBAR_GRAPH_TOP),
                sf::Vector2f(210.0f, SIDEBAR_GRAPH_HEIGHT)
            );
        }
        else if (graphMode == GraphMode::PhaseDiagram)
        {
            graphLabel.setString("Phase Diagram: Pressure vs Temp");
            window.draw(graphLabel);

            drawPhaseDiagram(
                window,
                phaseSamples,
                sf::Vector2f(temperature, pressure),
                sf::Vector2f(20.0f, SIDEBAR_GRAPH_TOP),
                sf::Vector2f(210.0f, SIDEBAR_GRAPH_HEIGHT)
            );
        }

        window.display();
    }

    return 0;
}


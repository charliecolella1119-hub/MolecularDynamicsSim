#include <SFML/Graphics.hpp>

#include <cmath>
#include <iostream>
#include <optional>
#include <random>
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
    Crystal
};

struct Particle
{
    sf::Vector2f position;
    sf::Vector2f velocity;
    sf::Vector2f force;
};

float lengthSquared(const sf::Vector2f& v)
{
    return v.x * v.x + v.y * v.y;
}

float speed(const Particle& p)
{
    return std::sqrt(lengthSquared(p.velocity));
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

float computeKineticEnergy(const std::vector<Particle>& particles)
{
    float kineticEnergy = 0.0f;

    for (const auto& p : particles)
    {
        kineticEnergy += 0.5f * MASS * lengthSquared(p.velocity);
    }

    return kineticEnergy;
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

            if (preset == Preset::Gas)
            {
                p.position = sf::Vector2f(
                    BOX_LEFT + spacingX * (x + 1),
                    BOX_TOP + spacingY * (y + 1)
                );
            }
            else if (preset == Preset::Liquid)
            {
                float clusterWidth = 320.0f;
                float clusterHeight = 260.0f;

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

    return "Crystal";
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
    graphLabel.setPosition(sf::Vector2f(25.0f, 405.0f));
    graphLabel.setString("Temperature History");

    sf::RectangleShape uiPanel;
    uiPanel.setPosition(sf::Vector2f(10.0f, 70.0f));
    uiPanel.setSize(sf::Vector2f(230.0f, 500.0f));
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
    std::vector<Particle> particles = createParticles(144, currentPreset);
    computeForces(particles);

    bool paused = false;

    std::vector<float> temperatureHistory;
    const size_t maxHistoryPoints = 180;

    float smoothedTemperature = 0.0f;
    bool hasSmoothedTemperature = false;

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
                    particles = createParticles(144, currentPreset);
                    computeForces(particles);
                    temperatureHistory.clear();
                    hasSmoothedTemperature = false;
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
                    currentPreset = Preset::Gas;
                    particles = createParticles(144, currentPreset);
                    computeForces(particles);
                    temperatureHistory.clear();
                    hasSmoothedTemperature = false;
                }

                if (keyPressed->code == sf::Keyboard::Key::Num2)
                {
                    currentPreset = Preset::Liquid;
                    particles = createParticles(144, currentPreset);
                    computeForces(particles);
                    temperatureHistory.clear();
                    hasSmoothedTemperature = false;
                }

                if (keyPressed->code == sf::Keyboard::Key::Num3)
                {
                    currentPreset = Preset::Crystal;
                    particles = createParticles(144, currentPreset);
                    computeForces(particles);
                    temperatureHistory.clear();
                    hasSmoothedTemperature = false;
                }
            }
        }

        if (!paused)
        {
            for (int i = 0; i < 4; i++)
            {
                velocityVerletStep(particles);
            }
        }

        float kineticEnergy = computeKineticEnergy(particles);
        float potentialEnergy = computePotentialEnergy(particles);
        float totalEnergy = kineticEnergy + potentialEnergy;
        float temperature =
            kineticEnergy / static_cast<float>(particles.size());

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

        infoText.setString(
            "Molecular Dynamics Simulator\n"
            "Preset: " + presetName(currentPreset) + "\n"
            "Temp: " + std::to_string(temperature) + "\n" +
            "Kinetic Energy: " + std::to_string(kineticEnergy) + "\n" +
            "Potential Energy: " + std::to_string(potentialEnergy) + "\n" +
            "Total Energy: " + std::to_string(totalEnergy) + "\n\n" +
            "Controls:\n"
            "1 = Gas\n"
            "2 = Liquid\n"
            "3 = Crystal\n"
            "Up = Heat\n"
            "Down = Cool\n"
            "Space = Pause\n"
            "R = Reset\n"
            "Esc = Quit"
        );

        window.setTitle(
            "MD Sim | Preset: " + presetName(currentPreset) +
            " | Temp: " + std::to_string(temperature)
        );

        window.clear(sf::Color(15, 17, 23));

        window.draw(box);

        for (const auto& p : particles)
        {
            particleShape.setPosition(p.position);
            particleShape.setFillColor(speedColor(speed(p)));
            window.draw(particleShape);
        }

        window.draw(uiPanel);
        window.draw(infoText);
        window.draw(graphLabel);

        drawTemperatureGraph(
            window,
            temperatureHistory,
            sf::Vector2f(20.0f, 430.0f),
            sf::Vector2f(210.0f, 120.0f)
        );

        window.display();
    }

    return 0;
}

